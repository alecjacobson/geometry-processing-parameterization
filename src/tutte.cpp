#include "tutte.h"
#include <Eigen/Sparse>
#include <vector>
#include <igl/boundary_loop.h>
#include <math.h>
#include <igl/min_quad_with_fixed.h>
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
typedef Eigen::Triplet<double> tuple;
void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  //space allocation
  int nV = V.rows();
  U.resize(nV, 2);

  //Extract indices of the boundary
  Eigen::VectorXi BndV; BndV.setZero();
  igl::boundary_loop(F, BndV);
  
  //UV coordinate around the circle
  //Add an offset term to align the texture
  int nBndV = BndV.rows();
  double offset = 22*M_PI/180;
  Eigen::MatrixXd UV; UV.setZero();
  UV.resize(nBndV, 2);
  for (int i = 0; i < nBndV; i++)
  {
    UV(i, 0) = sin(i*2*M_PI / nBndV + offset);
    UV(i, 1) = cos(i*2*M_PI / nBndV + offset);
  }

  //Compute L matrix
  int nF = F.rows();
  Eigen::SparseMatrix<double> L;
  L.resize(nV, nV); L.setZero();
  std::vector<tuple> tuple_list;
  for (int i = 0; i < nF; i++)
  {
    int iV0 = F(i, 0);
    int iV1 = F(i, 1);
    int iV2 = F(i, 2);
    double w0 = 1 / (V.row(iV0) - V.row(iV1)).norm();
    double w1 = 1 / (V.row(iV1) - V.row(iV2)).norm();
    double w2 = 1 / (V.row(iV2) - V.row(iV0)).norm();
    tuple_list.push_back(tuple(F(i, 0), F(i, 1), w0));
    tuple_list.push_back(tuple(F(i, 1), F(i, 2), w1));
    tuple_list.push_back(tuple(F(i, 2), F(i, 0), w2));
  }
  
  L.setFromTriplets(tuple_list.begin(), tuple_list.end());
  for (int i = 0; i < nV; i++)
    L.insert(i, i) = -L.row(i).sum();

  Eigen::SparseMatrix<double> A = -L;
  Eigen::VectorXd B = Eigen::VectorXd::Zero(nV, 1);
  Eigen::SparseMatrix<double> Aeq; Aeq.setZero();
  Eigen::MatrixXd Beq; Beq.setZero();
  igl::min_quad_with_fixed(A, B, BndV, UV, Aeq, Beq, false, U);
}