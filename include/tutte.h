#ifndef TUTTE_H
#define TUTTE_H
#include <Eigen/Core>
// Given a 3D mesh (`V`,`F`) with a disk topology (i.e., a manifold with single
// boundary), compute a 2D parameterization according to Tutte's mapping inside
// the unit disk. All boundary vertices should be mapped to the unit circle and
// interior vertices mapped inside the disk _without_ flips.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   U  #U by 2 list of mesh UV parameterization coordinates
//
void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U);
#endif
