#pragma once
#include<vector>
bool computeParametrization(
    const std::vector<std::array<double,3>>points,
    const std::vector<std::array<uint32_t,3>>& triangles,
    std::vector<std::array<double,3>>& stl_vertices_uv);   