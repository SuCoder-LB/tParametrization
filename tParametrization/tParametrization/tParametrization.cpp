#include "tParametrization.h"

#include<map>
#include<vector>
#include<array>
#include<algorithm>

#include"arrayGeometry.h"

#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/Dense>


using namespace ArrayGeometry;

bool SortEdgeConsecutive(std::vector<std::array<uint32_t,2>>& e,
    std::vector<std::vector<uint32_t> >& vs)
{
    vs.clear();
    if (e.empty()) return true;
    std::map<uint32_t, std::pair<uint32_t, uint32_t>> c;

    for (size_t i = 0; i < e.size(); i++) {
        uint32_t v0 = e[i][0];
        uint32_t v1 = e[i][1];

        auto it0 = c.find(v0), it1 = c.find(v1);
        if (it0 == c.end())
            c[v0] = std::make_pair(v1, (uint32_t)-1);
        else {
            if (it0->second.second == (uint32_t)-1) { it0->second.second = v1; }
            else {
                fprintf(stdout, "A list of edges has points that are adjacent to 3 edges\n");
                return false;
            }
        }
        if (it1 == c.end())
            c[v1] = std::make_pair(v0, (uint32_t)-1);
        else {
            if (it1->second.second == (uint32_t)-1) { it1->second.second = v0; }
            else {
                fprintf(stdout, "Wrong topology for a list of edges\n");
                fprintf(stdout, "Node %d is adjacent to more than 2 nodes %d %d\n",
                    v1, it1->second.first,it1->second.second);
                return false;
            }
        }
    }

    while (!c.empty()) {
        std::vector<uint32_t> v;
        uint32_t start = (uint32_t)-1;
        {
            auto it = c.begin();
            start = it->first;
            for (; it != c.end(); ++it) {
                if (it->second.second == (uint32_t)-1) {
                    start = it->first;
                    break;
                }
            }
        }

        auto its = c.find(start);

        uint32_t prev =
            (its->second.second == start) ? its->second.first : its->second.second;
        uint32_t current = start;

        do {
            if (c.size() == 0) {
                fprintf(stdout, "Wrong topology in a wire");
                return false;
            }
            v.push_back(current);
            auto it = c.find(current);
            if (it == c.end() || it->first == uint32_t(-1)) {
                fprintf(stdout, "Impossible to find %d", current);
                return false;
            }
            uint32_t v1 = it->second.first;
            uint32_t v2 = it->second.second;
            c.erase(it);
            uint32_t temp = current;
            if (v1 == prev)
                current = v2;
            else if (v2 == prev)
                current = v1;
            else {
                break;
            }
            prev = temp;
            if (current == start) { v.push_back(current); }
        } while (current != start && current != uint32_t(-1));
        if (v.size() > 2 && v[v.size() - 2] == v[v.size() - 1]) {
            v.erase(v.begin() + v.size() - 1);
        }
        vs.push_back(v);
    }
    return true;
}



bool computeParametrization(const std::vector<std::array<double, 3>> points, 
    const std::vector<std::array<uint32_t, 3>>& triangles,
    std::vector<std::array<double, 3>>& stl_vertices_uv)
{
    stl_vertices_uv.clear();


    if (triangles.empty()) return false;

    // get nodes and edges
    //std::map<MVertex*, int> nodeIndex;
    std::map<std::array<uint32_t, 2>, std::vector<std::array<uint32_t,3>>>edges;
    //std::map<MLine, std::vector<MTriangle*>, MLineLessThan> edges;
    for (std::size_t i = 0; i < triangles.size(); i++) {
        auto t = triangles[i];
        for (int j = 0; j < 3; j++) {
            
            std::array<uint32_t, 2>e = { t[j],t[(j + 1) % 3] };
            if (e[0] > e[1])std::swap(e[0], e[1]);
            edges[e].push_back(t);
        }
    }

    // compute edge loops
    std::vector<std::array<uint32_t,2>> es;
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        if (it->second.size() == 1) { // on boundary
            es.push_back(it->first);
        }
        else if (it->second.size() == 2) { // inside
        }
        else { // non-manifold: not supported
            fprintf(stdout, "Wrong topology of triangulation for parametrization: one edge is incident to %d triangles\n",
                it->second.size());
            return false;
        }
    }
    std::vector<std::vector<uint32_t> > vs;
    if (!SortEdgeConsecutive(es, vs)) {
        fprintf(stdout, "Wrong topology of boundary mesh for parametrization");
        return false;
    }
    if (vs.empty() || vs[0].size() < 2) {
        fprintf(stdout, "Invalid exterior boundary mesh for parametrization");
        return false;
    }

    //fprintf(stdout, "Parametrisation of surface with %lu triangles, %lu edges and %lu holes\n",
     //   triangles.size(), edges.size(), vs.size() - 1);

    // find longest loop and use it as the "exterior" loop
    int loop = 0;
    double longest = 0.;
    for (std::size_t i = 0; i < vs.size(); i++) {
        double l = 0.;
        for (std::size_t j = 1; j < vs[i].size(); j++) {
            l += length(points[vs[i][j]] - points[vs[i][j - 1]]);
        }
        if (l > longest) {
            longest = l;
            loop = i;
        }
    }

    // check orientation of the loop and reverse if necessary
    bool reverse = true;
    std::array<uint32_t, 2> ref = { vs[loop][0], vs[loop][1] };
    for (std::size_t i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        for (int j = 0; j < 3; j++) {
            if (t[j] == ref[0] && t[(j + 1) % 3] == ref[1]) {
                reverse = false;
                break;
            }
        }
        if (!reverse) break;
    }
    if (reverse) { std::reverse(vs[0].begin(), vs[0].end()); }

    std::vector<double> u(points.size(), 0.), v(points.size(), 0.);

    // boundary conditions
    std::vector<bool> bc(points.size(), false);
    double currentLength = 0;
    int index = vs[loop][0];
    bc[index] = true;
    u[index] = 1.;
    v[index] = 0.;
    for (std::size_t i = 1; i < vs[loop].size() - 1; i++) {
        currentLength += length(points[vs[loop][i]] - points[vs[loop][i - 1]]);
        double angle = 2 * M_PI * currentLength / longest;
        index = vs[loop][i];
        bc[index] = true;
        u[index] = cos(angle);
        v[index] = sin(angle);
    }

    Eigen::VectorXd X;
    Eigen::VectorXd B;
    Eigen::SparseMatrix<double> A;

    A.resize(points.size(), points.size());
    B.resize(points.size());
    X.resize(points.size());
    B.fill(0.);
    X.fill(0.);



    /*slow*/
    //{
    //    for (auto it = edges.begin(); it != edges.end(); ++it) {
    //        for (int ij = 0; ij < 2; ij++) {
    //            MVertex* v0 = it->first.getVertex(ij);
    //            int index0 = nodeIndex[v0];
    //            if (bc[index0]) continue; // boundary condition
    //            MVertex* v1 = it->first.getVertex(1 - ij);
    //            int index1 = nodeIndex[v1];
    //            MTriangle* tLeft = it->second[0];
    //            MVertex* vLeft = tLeft->getVertex(0);
    //            if (vLeft == v0 || vLeft == v1) vLeft = tLeft->getVertex(1);
    //            if (vLeft == v0 || vLeft == v1) vLeft = tLeft->getVertex(2);
    //            double e[3] = { v1->x() - v0->x(), v1->y() - v0->y(), v1->z() - v0->z() };
    //            double ne = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    //            double a[3] = { vLeft->x() - v0->x(), vLeft->y() - v0->y(),
    //                           vLeft->z() - v0->z() };
    //            double na = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    //            double thetaL =
    //                acos((a[0] * e[0] + a[1] * e[1] + a[2] * e[2]) / (na * ne));
    //            double thetaR = 0.;
    //            if (it->second.size() == 2) {
    //                MTriangle* tRight = it->second[1];
    //                MVertex* vRight = tRight->getVertex(0);
    //                if (vRight == v0 || vRight == v1) vRight = tRight->getVertex(1);
    //                if (vRight == v0 || vRight == v1) vRight = tRight->getVertex(2);
    //                double b[3] = { vRight->x() - v0->x(), vRight->y() - v0->y(),
    //                               vRight->z() - v0->z() };
    //                double nb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    //                thetaR = acos((b[0] * e[0] + b[1] * e[1] + b[2] * e[2]) / (nb * ne));
    //            }
    //            double c = (tan(.5 * thetaL) + tan(.5 * thetaR)) / ne;
    //            lsys->addToMatrix(index0, index1, -c);
    //            lsys->addToMatrix(index0, index0, c);
    //        }
    //    }
    //    for (std::size_t i = 0; i < vs[loop].size() - 1; i++) {
    //        int row = nodeIndex[vs[loop][i]];
    //        lsys->addToMatrix(row, row, 1);
    //    }
    //}


    /*faster*/
    {
        std::vector<Eigen::Triplet<double, size_t> > triplets;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            for (int ij = 0; ij < 2; ij++) {
                uint32_t v0 = it->first[ij];
                if (bc[v0]) continue; // boundary condition
                uint32_t v1 = it->first[1-ij];
               
                auto tLeft = it->second[0];
                uint32_t vLeft = tLeft[0];
                if (vLeft == v0 || vLeft == v1) vLeft = tLeft[1];
                if (vLeft == v0 || vLeft == v1) vLeft = tLeft[2];
                vec3 e = points[v1] - points[v0];
         
                double ne = length(e);
                vec3 a = points[vLeft] - points[v0];
    
                double thetaL = angleVectors(e, a);
                
                double thetaR = 0.;
                if (it->second.size() == 2) {
                    auto& tRight = it->second[1];
                    uint32_t vRight = tRight[0];
                    if (vRight == v0 || vRight == v1) vRight = tRight[1];
                    if (vRight == v0 || vRight == v1) vRight = tRight[2];
                    vec3 b = points[vRight] - points[v0];
                    thetaR = angleVectors(e, b);
            
                }
                double c = (tan(.5 * thetaL) + tan(.5 * thetaR)) / ne;
                triplets.push_back({ (size_t)v0,(size_t)v1,-c });
                triplets.push_back({ (size_t)v0,(size_t)v0,c });
            }
        }
        for (std::size_t i = 0; i < vs[loop].size() - 1; i++) {
            int row = vs[loop][i];
            triplets.push_back({ (size_t)row,(size_t)row,1 });
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
    }


    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);


    B.fill(0.);
    for (std::size_t i = 0; i < vs[loop].size() - 1; i++) {
        int row = vs[loop][i];
        B[row] += u[row];
    }

    X = solver.solve(B);


    for (std::size_t i = 0; i < points.size(); i++) {
        u[i] = X[i];
    }

    B.fill(0.);
    for (std::size_t i = 0; i < vs[loop].size() - 1; i++) {
        int row = vs[loop][i];
        B[row] += v[row];
    }
    X = solver.solve(B);

    for (std::size_t i = 0; i < points.size(); i++) {
        v[i] = X[i];
    }

    stl_vertices_uv.resize(points.size());

    for (std::size_t i = 0; i < points.size(); i++) {
        stl_vertices_uv[i] = vec3{ u[i], v[i],0 };
    }

    return true;
}
