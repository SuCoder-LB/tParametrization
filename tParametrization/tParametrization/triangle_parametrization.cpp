#include "triangle_parametrization.h"

#include <numeric>

#include "basic_types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

bool SortEdgeConsecutive(const std::vector<ID2>& e,
  std::vector<std::vector<ID> >& vs) {
  vs.clear();
  if (e.empty()) return true;
  std::map<ID, std::pair<ID, ID>> c;

  for (auto i : e) {
    ID v0 = i[0];
    ID v1 = i[1];
    auto it0 = c.find(v0), it1 = c.find(v1);
    if (it0 == c.end())
      c[v0] = std::make_pair(v1, NO_ID);
    else {
      if (it0->second.second == NO_ID) { it0->second.second = v1; }
      else {
        debug(__FILE__,
          __LINE__,
          "A list of edges has points that are adjacent to 3 edges");
        return false;
      }
    }
    if (it1 == c.end())
      c[v1] = std::make_pair(v0, NO_ID);
    else {
      if (it1->second.second == NO_ID) { it1->second.second = v0; }
      else {
        debug(__FILE__, __LINE__, "Wrong topology for a list of edges");
        debug(__FILE__,
          __LINE__,
          "Node %d is adjacent to more than 2 nodes %d %d",
          v1,
          it1->second.first,
          it1->second.second);
        return false;
      }
    }
  }

  while (!c.empty()) {
    std::vector<ID> v;
    ID start = NO_ID;
    {
      auto it = c.begin();
      start = it->first;
      for (; it != c.end(); ++it) {
        if (it->second.second == NO_ID) {
          start = it->first;
          break;
        }
      }
    }

    auto its = c.find(start);

    ID prev =
      (its->second.second == start) ? its->second.first : its->second.second;
    ID current = start;

    do {
      if (c.empty()) {
        warn(__FILE__, __LINE__, "Wrong topology in a wire");
        return false;
      }
      v.push_back(current);
      auto it = c.find(current);
      if (it == c.end() || it->first == NO_ID) {
        error(__FILE__, __LINE__, "Impossible to find %d", current);
        return false;
      }
      ID v1 = it->second.first;
      ID v2 = it->second.second;
      c.erase(it);
      ID temp = current;
      if (v1 == prev)
        current = v2;
      else if (v2 == prev)
        current = v1;
      else {
        break;
      }
      prev = temp;
      if (current == start) { v.push_back(current); }
    } while (current != start && current != NO_ID);
    if (v.size() > 2 && v[v.size() - 2] == v[v.size() - 1]) {
      v.erase(v.begin() + (int)v.size() - 1);
    }
    vs.push_back(v);
  }
  return true;
}

template<class T>
bool buildBoundaries(const std::vector<T>& elements,
  std::vector<std::vector<uint32_t> >& loops) {
  loops.clear();

  std::vector<std::array<uint32_t, 2>> eds, v_eds;
  for (auto e : elements) {
    for (int j = 0;j < (int)e.size();++j) {
      eds.push_back({ e[j], e[(j + 1) % e.size()] });
      std::sort(eds.back().begin(), eds.back().end());
    }
  }

  std::sort(eds.begin(), eds.end());

  for (size_t i = 0; i < eds.size(); i++) {
    if (i != eds.size() - 1 && eds[i] == eds[i + 1])i++;
    else v_eds.push_back(eds[i]);
  }

  if (v_eds.empty()) return true; /* No boundary, e.g. sphere */

  std::vector<std::vector<ID> > v_sorted;
  bool oks = SortEdgeConsecutive(v_eds, v_sorted);
  if (!oks || v_sorted.empty()) {
    return false;
  }

  /* Reverse vertices if necessary, to keep coherent with elements orientation */
  for (auto l : v_sorted) {
    std::vector<ID>& loop = l;
    /* Find a MEdge on the loop */
    auto a = (ID)-1;
    auto b = (ID)-1;
    for (auto e : v_eds) {
      ID v1 = e[0];
      ID v2 = e[1];
      auto it = std::find(loop.begin(), loop.end(), v1);
      if (it != loop.end()) {
        /* v1 from the MEdge found on the loop */
        a = v1;
        b = v2;
        break;
      }
    }
    if (a == (uint32_t)-1 || b == (uint32_t)-1) {
      fprintf(stdout,
        "buildBoundaries(): vertex not found in Sorted vertices, weird\n");
      return false;
    }

    auto it = std::find(loop.begin(), loop.end(), a);
    size_t i = it - loop.begin();
    size_t i_next = (i + 1) % loop.size();
    size_t i_prev = (i - 1 + loop.size()) % loop.size();
    if (loop[i_next] == b) {
      // good ordering
    }
    else if (loop[i_prev] == b) { // apply reverse
      std::reverse(loop.begin(), loop.end());
    }
    else {
      fprintf(stdout,
        "buildBoundaries(): second vertex not found in adjacent Sorted vertices, weird\n");
      return false;
    }
  }

  loops = v_sorted;
  return true;
}

bool computeParametrization(const std::vector<std::array<double, 3>>& points,
  const std::vector<std::array<uint32_t,
  3>> &triangles,
  std::vector<std::array<double,
  3>> *stl_vertices_uv) {
  (*stl_vertices_uv).clear();
  (*stl_vertices_uv).resize(points.size());
  if (triangles.empty()) return false;

  std::vector<std::vector<uint32_t> > vs;
  if (!buildBoundaries(triangles, vs))return false;

  // find the longest loop and use it as the "exterior" loop
  int loop = 0;
  double longest = 0.;
  for (int i = 0; i < (int)vs.size(); i++) {
    double l = 0.;
    for (int j = 1; j < (int)vs[i].size(); j++) {
      l += Length(points[vs[i][j]] - points[vs[i][j - 1]]);
    }
    if (l > longest) {
      longest = l;
      loop = i;
    }
  }

  std::vector<double> b_u(points.size(), 0.), b_v(points.size(), 0.);

  // boundary conditions
  std::vector<bool> bc(points.size(), false);
  double currentLength = 0;

  for (int i = 0; i < (int)vs[loop].size(); i++) {
    if (i >= 1)
      currentLength += Length(points[vs[loop][i]] - points[vs[loop][i - 1]]);
    double angle = 2 * M_PI * currentLength / longest;
    int index = (int)vs[loop][i];
    bc[index] = true;
    b_u[index] = cos(angle);
    b_v[index] = sin(angle);
  }

  Eigen::VectorXd x_;
  Eigen::VectorXd b_;
  Eigen::SparseMatrix<double> a_;

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_lu;


  a_.resize((int)points.size(), (int)points.size());
  b_.resize((int)points.size());
  x_.resize((int)points.size());
  b_.fill(0.);
  x_.fill(0.);

  std::vector<std::vector<ID>> v2v;
  v2v.resize(points.size());
  for (auto& f : triangles) {
    for (size_t le = 0; le < f.size(); ++le) {
      ID vst[2] = { f[le], f[(le + 1) % f.size()] };
      if (vst[0] < vst[1]) {
        v2v[vst[0]].push_back(vst[1]);
        v2v[vst[1]].push_back(vst[0]);
      }
    }
  }

  std::vector<std::vector<size_t>> column(points.size());
  std::vector<std::vector<double>> value(points.size());
  for (size_t p = 0; p < points.size(); ++p) {
    if (!bc[p]) {
      column[p].push_back(p);
      value[p].push_back(1);

      if (v2v[p].empty()) continue;
      auto sum = double(v2v[p].size());
      for (size_t v2 : v2v[p]) {
        column[p].push_back(v2);
        value[p].push_back(-1. / sum);
      }
    }
    else { /* fixed value */
      column[p].push_back((int)p);
      value[p].push_back(1);
    }
  }


  // Msg::Debug("Eigen call: add coefficients");
  std::vector<Eigen::Triplet<double, size_t> > triplets;
  triplets.reserve(value.size());

  for (size_t i = 0; i < column.size(); ++i) {
    for (size_t k = 0; k < column[i].size(); ++k) {
      size_t j = column[i][k];
      double val = value[i][k];
      triplets.emplace_back(Eigen::Triplet<double, size_t>{i, j, val});
    }
  }
  a_.setFromTriplets(triplets.begin(), triplets.end());
  
  bool solveOk = true;
  { /* Try SparseLU */
    solver_lu.analyzePattern(a_);
    DBG("preprocess");
    solver_lu.factorize(a_);
    DBG("factorize");

    for (int i = 0; i < (int)b_u.size(); ++i) {
      if (b_u[i] != 0.) {
        b_[i] = b_u[i];
      }
    }

    x_ = solver_lu.solve(b_);
    
    b_u.clear();
    b_u.resize(x_.size());
    for (int i = 0; i < x_.size(); ++i) b_u[i] = x_[i];


    for (int i = 0; i < (int)b_v.size(); ++i) {
      if (b_v[i] != 0.) {
        b_[i] = b_v[i];
      }
    }

    x_ = solver_lu.solve(b_);

    b_v.clear();
    b_v.resize(x_.size());
    for (int i = 0; i < x_.size(); ++i) b_v[i] = x_[i];
  }

  for (size_t v = 0; v < points.size(); ++v) {
    (*stl_vertices_uv)[v][0] = b_u[v];
    (*stl_vertices_uv)[v][1] = b_v[v];
  }
  if (!solveOk) {
    error(__FILE__, __LINE__, "failed to solve linear system to solve uv");
    return false;
  }

  return true;
}

void InitUV(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<Vec3>* stl_vertices_uv) {
  std::vector<Eigen::Vector2d> uvs;
  computeParametrization(points, triangles, stl_vertices_uv);
}

std::vector<std::vector<Vec3>> CalEdgeVectors(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles) {

  std::vector<std::vector<Vec3>> vectors;
  vectors.resize(triangles.size());

  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {

    Vec3 arr[3] = { points[triangles[f_it][0]], points[triangles[f_it][1]],
                   points[triangles[f_it][2]] };

    double a = Length(arr[0] - arr[1]);
    double b = Length(arr[1] - arr[2]);
    double c = Length(arr[2] - arr[0]);
    double angle = std::acos((a * a + c * c - b * b) / (2 * a * c));
    Vec3 UPosition1{ 0, 0, 0 };
    Vec3 UPosition2{ a, 0, 0 };
    Vec3 UPosition3{ c * cos(angle), c * sin(angle), 0 };

    vectors[f_it].push_back(UPosition3 - UPosition2);
    vectors[f_it].push_back(UPosition1 - UPosition3);
    vectors[f_it].push_back(UPosition2 - UPosition1);
  }
  return vectors;
}

double** CalCots(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles) {
  double** C;
  C = new double* [triangles.size()];
  for (int i = 0; i < (int)triangles.size(); ++i) {
    C[i] = new double[3];
  }

  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {

    for (int i = 0; i < 3; ++i) {
      Vec3 prep =
        points[triangles[f_it][(i + 2) % 3]] - points[triangles[f_it][i]];
      Vec3 next_p =
        points[triangles[f_it][(i + 1) % 3]] - points[triangles[f_it][i]];
      Eigen::Vector3d pre = { prep[0], prep[1], prep[2] };
      Eigen::Vector3d next = { next_p[0], next_p[1], next_p[2] };
      pre.normalize();
      next.normalize();
      double angle = acos(pre.dot(next));
      C[f_it][i] = 1.0 / std::tan(angle);
    }
  }
  return C;
}
double CalRigidEnergy(const std::vector<ID3>& triangles,
  std::vector<std::vector<Vec3>>& EV,
  std::vector<Vec3>& vt,
  double** C,
  std::vector<Eigen::Matrix2d> R) {
  double E = 0.0;
  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {
    for (int i = 0; i < 3; ++i) {
      int n = (i + 1) % 3;
      int p = (i + 2) % 3;
      Vec3 Uij = vt[triangles[f_it][p]] - vt[triangles[f_it][n]];
      Vec3 Eij = EV[f_it][i];
      Vec3 REij = { R[f_it](0, 0) * Eij[0] + R[f_it](0, 1) * Eij[1],
                   R[f_it](1, 0) * Eij[0] + R[f_it](1, 1) * Eij[1] };
      E += C[f_it][i] * (std::pow(Length(Uij - REij), 2));
      //E += std::pow(C[f_it][i] * (Uij - R[f_it] * Eij).norm(), 2);
    }
  }
  return E;
}
Eigen::MatrixXd ComputeTransformMatrix(
  const std::vector<ID3>& triangles,
  const std::vector<Vec3>& vt,
  const std::vector<std::vector<Vec3>>& EV,
  double** C,
  int face_id) {
  std::vector<Eigen::Vector2d> u;

  for (int i = 0; i < 3; ++i)
    u
    .emplace_back(Eigen::Vector2d(vt[triangles[face_id][i]][0],
      vt[triangles[face_id][i]][1]));

  Eigen::MatrixXd uu, xx, cc;
  uu.resize(3, 3);
  xx.resize(3, 3);
  cc.resize(3, 3);
  for (int i = 0; i < 3; ++i) {
    uu(i, 0) = 0;
    uu(i, 1) = 0;
    uu(i, 2) = 0;
    xx(i, 0) = 0;
    xx(i, 1) = 0;
    xx(i, 2) = 0;
    cc(i, 0) = 0;
    cc(i, 1) = 0;
    cc(i, 2) = 0;
  }

  for (int i = 0; i < 3; i++) {
    int n = (i + 1) % 3;
    int p = (i + 2) % 3;
    uu(i, 0) = u[p][0] - u[n][0];
    uu(i, 1) = u[p][1] - u[n][1];
    xx(i, 0) = EV[face_id][i][0];
    xx(i, 1) = EV[face_id][i][1];
    cc(i, i) = C[face_id][i];
  }
  Eigen::MatrixXd CovMat = xx.transpose() * cc * uu;


  /*Eigen::MatrixXd uu1, xx1, cc1;
  uu1.resize(3, 3); xx1.resize(3, 3);

  for (int i = 0; i < 2; ++i) {
      uu1(i, 0) = 0; uu1(i, 1) = 0; uu1(i, 2) = 0;
      xx1(i, 0) = 0; xx1(i, 1) = 0; xx1(i, 2) = 0;
  }
  uu1(2, 0) = 1; uu1(2, 1) = 1; uu1(2, 2) = 1;
  xx1(2, 0) = 1; xx1(2, 1) = 1; xx1(2, 2) = 1;
  for (int i = 0; i < 2; i++) {
      int n = (i + 1) % 3;
      int p = (i + 2) % 3;
      uu1(0, i) = u[p][0] - u[n][0];
      uu1(1, i) = u[p][1] - u[n][1];
      xx1(0, i) = EV[face_id][i][0];
      xx1(1, i) = EV[face_id][i][1];
  }
  Eigen::MatrixXd CovMat1 = uu1 * xx1.inverse();*/

  return CovMat;
  //return CovMat1.transpose();
}
Eigen::MatrixXd SVDFactorize(const Eigen::MatrixXd& CovMat) {

  Eigen::BDCSVD<Eigen::MatrixXd>
    svd(CovMat, Eigen::ComputeThinU | Eigen::ComputeThinV);

  Eigen::MatrixXd matrixU;
  matrixU.resize(3, 3);
  matrixU = svd.matrixU();
  Eigen::MatrixXd rot = svd.matrixV() * (matrixU.transpose());
  if (rot.determinant() < 0) {
    if (svd.singularValues()[0] < svd.singularValues()[1]) {
      for (int i = 0; i < 2; ++i) {
        matrixU(i, 0) = -1 * svd.matrixU()(i, 0);
      }
    }
    else {
      for (int i = 0; i < 2; ++i) {
        matrixU(i, 1) = -1 * svd.matrixU()(i, 1);
      }
    }
    rot = svd.matrixV() * (matrixU.transpose());

  }
  return rot;
}
using namespace std;

std::vector<Eigen::Matrix2d> ARAPLocal(const std::vector<ID3>& triangles,
  const std::vector<Vec3>& vt,
  std::vector<std::vector<Vec3>>& EV,
  double** C,
  std::vector<Eigen::Matrix2d>& R) {
  R.clear();
  R.resize(triangles.size());

  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {
    Eigen::MatrixXd CovMat = ComputeTransformMatrix(triangles,
      vt,
      EV,
      C,
      f_it);
    Eigen::MatrixXd rot = SVDFactorize(CovMat);
    R[f_it] << rot(0, 0), rot(0, 1), rot(1, 0), rot(1, 1);
  }
  return R;
}

bool ARAPGlobal(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<std::vector<Vec3>>& EV,
  Eigen::SparseMatrix<double>& L,
  double** C,
  std::vector<Eigen::Matrix2d>& R,
  Eigen::SparseLU<Eigen::SparseMatrix<double> >& Solver,
  std::vector<Vec3>* vt) {

  int l = (int)points.size();
  Eigen::VectorXd bx;
  bx.resize(l);
  Eigen::VectorXd by;
  by.resize(l);

  for (int i = 0; i < l; ++i) {
    bx[i] = 0;
    by[i] = 0;
  }

  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {
    for (int i = 0; i < 3; ++i) {
      int n = (i + 1) % 3;
      int p = (i + 2) % 3;
      Vec3 Eij = EV[f_it][i];
      double x = R[f_it](0, 0) * Eij[0] * C[f_it][i]
        + R[f_it](0, 1) * Eij[1] * C[f_it][i];
      double y = R[f_it](1, 0) * Eij[0] * C[f_it][i]
        + R[f_it](1, 1) * Eij[1] * C[f_it][i];
      bx[triangles[f_it][n]] -= x;
      bx[triangles[f_it][p]] += x;
      by[triangles[f_it][n]] -= y;
      by[triangles[f_it][p]] += y;
    }
  }
  //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> Solver1, Solver2;
  //Solver1.compute(L); Solver2.compute(L);
  fprintf(stdout, "Eigen call: solve linear system\n");
  Eigen::VectorXd Ux = Solver.solve(L * bx);
  if (Solver.info() != Eigen::ComputationInfo::Success) {
    fprintf(stdout, "Eigen: failed to solve linear system with SparseLU\n");
    return false;
  }
  fprintf(stdout, "Eigen call: solve linear system\n");
  Eigen::VectorXd Uy = Solver.solve(L * by);
  if (Solver.info() != Eigen::ComputationInfo::Success) {
    fprintf(stdout, "Eigen: failed to solve linear system with SparseLU\n");
    return false;
  }

  //Eigen::VectorXd Ux = Solver.solve(bx);
  //Eigen::VectorXd Uy = Solver.solve(by);

  for (int i = 0; i < l; ++i) {
    (*vt)[i][0] = Ux[i];
    (*vt)[i][1] = Uy[i];
  }
  return true;
}

void Laplacian(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<std::vector<size_t>>* column,
  std::vector<std::vector<double>>* value) {

  for (auto t : triangles) {

    Vec3 v1 = points[t[0]];
    Vec3 v2 = points[t[1]];
    Vec3 v3 = points[t[2]];

    double cot1 = Dot(v2 - v1, v3 - v1) / Length(Cross(v2 - v1, v3 - v1));
    double cot2 = Dot(v3 - v2, v1 - v2) / Length(Cross(v3 - v2, v1 - v2));
    double cot3 = Dot(v1 - v3, v2 - v3) / Length(Cross(v1 - v3, v2 - v3));
    (*column)[t[0]].push_back(t[1]);
    (*value)[t[0]].push_back(cot3);
    (*column)[t[1]].push_back(t[0]);
    (*value)[t[1]].push_back(cot3);

    (*column)[t[1]].push_back(t[2]);
    (*value)[t[1]].push_back(cot1);
    (*column)[t[2]].push_back(t[1]);
    (*value)[t[2]].push_back(cot1);

    (*column)[t[2]].push_back(t[0]);
    (*value)[t[2]].push_back(cot2);
    (*column)[t[0]].push_back(t[2]);
    (*value)[t[0]].push_back(cot2);

  }

  for (int i = 0; i < (int)points.size(); ++i) {
    (*column)[i].push_back(i);
    double sum = std::accumulate((*value)[i].begin(), (*value)[i].end(), 0.0);
    (*value)[i].push_back(sum * -1);
  }

}

using namespace std;

bool ComputeAsRigidAsPossibleParametrization(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<Vec3>* stl_vertices_uv) {
  InitUV(points, triangles, stl_vertices_uv);
  DBG("init success");

  std::vector<std::vector<Vec3>> EV = CalEdgeVectors(points, triangles);
  double** C = CalCots(points, triangles);
  DBG("cal cots success");
  Eigen::SparseMatrix<double> L((int)points.size(), (int)points.size());
  std::vector<std::vector<size_t>> columns(points.size());
  std::vector<std::vector<double>> values(points.size());
  Laplacian(points, triangles, &columns, &values);

  DBG("laplacian");
  std::vector<Eigen::Triplet<double, size_t> > triplets;
  triplets.reserve(values.size());

  for (size_t i = 0; i < columns.size(); ++i) {
    for (size_t k = 0; k < columns[i].size(); ++k) {
      size_t j = columns[i][k];
      double val = values[i][k];
      triplets.emplace_back(Eigen::Triplet<double, size_t>{i, j, val});
    }
  }
  L.setFromTriplets(triplets.begin(), triplets.end());

  DBG("laplacian");

  double E = -1;
  double energy_pre = 0;
  int iterations = 0;
  Eigen::SparseMatrix<double> L2 = L * L;

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  debug(__FILE__,
    __LINE__,
    "Eigen call: analyse sparse matrix sparsity pattern\n");
  solver.analyzePattern(L2);

  debug(__FILE__, __LINE__, "Eigen call: factorize sparse matrix\n");
  solver.factorize(L2);

  if (solver.info() != Eigen::ComputationInfo::Success) {
    error(__FILE__,
      __LINE__,
      "Eigen: failed to solve linear system with SparseLU\n");
    return -1;
  }
  std::vector<Eigen::Matrix2d> R;

  while (std::abs(energy_pre - E) > 1e-6 && iterations < 30) {
    ++iterations;
    energy_pre = E;

    ARAPLocal(triangles, *stl_vertices_uv, EV, C, R);
    ARAPGlobal(points, triangles, EV, L, C, R, solver, stl_vertices_uv);

    E = CalRigidEnergy(triangles, EV, *stl_vertices_uv, C, R);
    std::cout << "the iterations is   :   " << iterations << std::endl;
    std::cout << "the energy is   :   " << E << std::endl << std::endl;

  }

  return true;
}

int ComputerLeastSquareConformalMaps(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<Vec3>* stl_vertices_uv,
  int fix1,
  int fix2) {
  stl_vertices_uv->clear();
  stl_vertices_uv->resize(points.size(), { 0, 0, 0 });
  /*初始化边长向量，用于后续构建矩阵*/
  std::vector<std::vector<Vec3> > vectors;
  vectors.resize(triangles.size());
  std::vector<Vec3> vertexs(3);
  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {
    for (int id = 0; id < 3; ++id)vertexs[id] = points[triangles[f_it][id]];
    double a = Length(vertexs[0] - vertexs[1]);
    double b = Length(vertexs[1] - vertexs[2]);
    double c = Length(vertexs[2] - vertexs[0]);
    double angle = std::acos((a * a + c * c - b * b) / (2 * a * c));
    Vec3 UPosition1{ 0, 0, 0 };
    Vec3 UPosition2{ a, 0, 0 };
    Vec3 UPosition3{ c * cos(angle), c * sin(angle), 0 };
    vectors[f_it].push_back(UPosition3 - UPosition2);
    vectors[f_it].push_back(UPosition1 - UPosition3);
    vectors[f_it].push_back(UPosition2 - UPosition1);
  }
  /*定义两个数组保证正向反向映射序号*/
  std::vector<int> VertexMapping(points.size());
  std::vector<int> antiVertexMapping(points.size());
  VertexMapping[fix1] = (int)points.size() - 2;
  VertexMapping[fix2] = (int)points.size() - 1;
  antiVertexMapping[points.size() - 2] = fix1;
  antiVertexMapping[points.size() - 1] = fix2;
  int itrcount = 0;
  for (int v_it = 0; v_it < (int)points.size(); ++v_it) {
    if (v_it == fix1 || v_it == fix2) continue;
    VertexMapping[v_it] = itrcount;
    antiVertexMapping[itrcount] = v_it;
    ++itrcount;
  }
  DBG("here");
  //构建A,B矩阵
  std::vector<Eigen::Triplet<double>> tripletList1;
  std::vector<Eigen::Triplet<double>> tripletList2;
  for (int f_it = 0; f_it < (int)triangles.size(); ++f_it) {
    double k = 1.0 / std::sqrt(
      (vectors[f_it][0][1]) * (vectors[f_it][2][0]) / 2.0);

    for (int fv_it = 0; fv_it < 3; ++fv_it) {
      if ((int)triangles[f_it][fv_it] == fix1
        || (int)triangles[f_it][fv_it] == fix2) {
        tripletList2.emplace_back(Eigen::Triplet<double>(f_it,
          VertexMapping[triangles[f_it][fv_it]]
          - (int)points
          .size() + 2,
          vectors[f_it][fv_it][0]
          * k));
        tripletList2.emplace_back(Eigen::Triplet<double>(f_it,
          VertexMapping[triangles[f_it][fv_it]]
          - (int)points
          .size() + 4,
          -vectors[f_it][fv_it][1]
          * k));
        tripletList2
          .emplace_back(Eigen::Triplet<double>(f_it + (int)triangles.size(),
            VertexMapping[triangles[f_it][fv_it]]
            - (int)points.size() + 2,
            vectors[f_it][fv_it][1] * k));
        tripletList2
          .emplace_back(Eigen::Triplet<double>(f_it + (int)triangles.size(),
            VertexMapping[triangles[f_it][fv_it]]
            - (int)points.size() + 4,
            vectors[f_it][fv_it][0] * k));
      }
      else {
        tripletList1.emplace_back(Eigen::Triplet<double>(f_it,
          VertexMapping[triangles[f_it][fv_it]],
          vectors[f_it][fv_it][0]
          * k));
        tripletList1.emplace_back(Eigen::Triplet<double>(f_it,
          VertexMapping[triangles[f_it][fv_it]]
          + (int)points
          .size() - 2,
          -vectors[f_it][fv_it][1]
          * k));
        tripletList1
          .emplace_back(Eigen::Triplet<double>(f_it + (int)triangles.size(),
            VertexMapping[triangles[f_it][fv_it]],
            vectors[f_it][fv_it][1] * k));
        tripletList1
          .emplace_back(Eigen::Triplet<double>(f_it + (int)triangles.size(),
            VertexMapping[triangles[f_it][fv_it]]
            + (int)points.size() - 2,
            vectors[f_it][fv_it][0] * k));
      }
    }
  }

  Eigen::SparseMatrix<double>
    A(2 * (int)triangles.size(), 2 * (int)points.size() - 4);
  Eigen::SparseMatrix<double> B(2 * (int)triangles.size(), 4);
  A.setFromTriplets(tripletList1.begin(), tripletList1.end());
  B.setFromTriplets(tripletList2.begin(), tripletList2.end());

  //  cout<<"["<<endl;
  //  for(int i=0;i<A.rows();++i){
  //    for(int j=0;j<A.cols();++j){
  //      cout<<A.coeff(i,j)<<" ";
  //    }
  //    cout<<";"<<endl;
  //  }
  //  cout<<"]"<<endl;
    //DBG(A,B);

    /*设置固定点的坐标*/
  Eigen::Vector4d U;
  U[0] = 0;
  U[1] = 1;
  U[2] = 0;
  U[3] = 0;
  /*计算右侧矩阵*/
  Eigen::VectorXd BM = -B * U;
  DBG("here");
  auto A2 = A.transpose() * A;

  DBG("A2 is ok");

  auto ABM = A.transpose() * BM;

  DBG("ABM is ok");
  /*求解线性方程组*/
  Eigen::SparseLU<Eigen::SparseMatrix<double>> Solver;
  Solver.analyzePattern(A2);
  Solver.factorize(A2);

  DBG("here");
  DBG(A.rows(), A.cols(), BM.rows());

  Eigen::VectorXd newP(2 * (int)points.size() - 4);
  newP = Solver.solve(ABM);

  DBG("here");
  /*把得到的新坐标映射回原来的顺序*/
  (*stl_vertices_uv)[fix1][0] = U[0];
  (*stl_vertices_uv)[fix1][1] = U[2];
  (*stl_vertices_uv)[fix2][0] = U[1];
  (*stl_vertices_uv)[fix2][1] = U[3];
  DBG("here");
  DBG(newP.rows());
  //DBG(newP);

  for (int i = 0; i < (int)VertexMapping.size() - 2; ++i) {
    int index = antiVertexMapping[i];
    //DBG(i, i + points.size() - 2, index);
    (*stl_vertices_uv)[index][0] = newP(i);
    (*stl_vertices_uv)[index][1] = newP(i + (int)points.size() - 2);
  }
  return 1;
}

int ComputerLeastSquareConformalMaps(const std::vector<Vec3>& points,
  const std::vector<ID3>& triangles,
  std::vector<Vec3>* stl_vertices_uv) {
  return ComputerLeastSquareConformalMaps(points, triangles, stl_vertices_uv, 0, (int)points.size() / 2);
}

