#pragma once

#include <vector>
#include <array>
#include <cstdint>

bool computeParametrization(const std::vector<std::array<double, 3>>& points,
  const std::vector<std::array<uint32_t,
  3>> &triangles,
  std::vector<std::array<double,
  3>> *stl_vertices_uv);

bool ComputeAsRigidAsPossibleParametrization(const std::vector<std::array<double,
  3>> &verts,
  const std::vector<std::array<
  uint32_t,
  3>> &triangles,
  std::vector<std::array<double,
  3>> *stl_vertices_uv);


int ComputerLeastSquareConformalMaps(const std::vector<std::array<double, 3>>& verts,
  const std::vector<std::array<uint32_t, 3>>& triangles,
  std::vector<std::array<double, 3>>* stl_vertices_uv);