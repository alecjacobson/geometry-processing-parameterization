#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/edges.h>
#include <igl/cotmatrix.h>
#include <iostream>
using namespace std;

// =====
// RIP
// =====
// void compStiffMat(
// 	const Eigen::MatrixXd & V,
// 	const Eigen::MatrixXi & F,
// 	Eigen::SparseMatrix<double> & L){

// 	// get edge list
// 	Eigen::MatrixXi E;
// 	igl::edges(F, E);

// 	// assemble tripletlist
// 	vector< Eigen::Triplet<double> > tripletList;
// 	tripletList.reserve(E.rows() * 3); // "E.rows() * 3" is the estimation_of_entries
// 	Eigen::VectorXd diagVals; // for storing diagonal values
// 	diagVals.setZero(V.rows());
// 	for (int ii = 0; ii < E.rows(); ii++){ // non-diagonal terms
// 		int VaIdx = E(ii,0);
// 		int VbIdx = E(ii,1);
// 		double weight = 1 / (V.row(VaIdx) - V.row(VbIdx)).norm();
// 		tripletList.emplace_back(VaIdx, VbIdx, weight);
// 		tripletList.emplace_back(VbIdx, VaIdx, weight);
// 		diagVals(VaIdx) -= weight;
// 		diagVals(VbIdx) -= weight;
// 	}
// 	for (int ii = 0; ii < V.rows(); ii++){ // diagonal terms
// 		tripletList.emplace_back(ii,ii,diagVals(ii));
// 	}

// 	// assemble stiffness matrix
// 	L.resize(V.rows(),V.rows());
// 	L.setFromTriplets(tripletList.begin(),tripletList.end());
// }

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = V.leftCols(2);

	// My code
	// extract boundary vertex indices
	vector<vector<int> > bdLoops;
	igl::boundary_loop(F,bdLoops);

	// only pin down the longest bd loop
	Eigen::VectorXi longestBdLoop;
	longestBdLoop.resize(bdLoops[0].size()); 
	for (int ii = 0; ii < bdLoops[0].size(); ii++){ 
		longestBdLoop(ii) = bdLoops[0][ii];
	}

	// pin boundary to circle
	Eigen::MatrixXd VCir;
	igl::map_vertices_to_circle(V,longestBdLoop,VCir);

	// minimize quaratic energy with constraints on the solution
	Eigen::SparseMatrix<double> L; // stiffness matrix
	// compStiffMat(V,F,L);
	igl::cotmatrix(V,F,L); // we can use cotanLaplace

	Eigen::SparseMatrix<double> Aeq; // empty matrix
	Eigen::VectorXd Beq; // empty matrix
	Eigen::VectorXd B;
	B.setZero(V.rows());

	igl::min_quad_with_fixed(L, B, longestBdLoop, VCir, Aeq, Beq, false, U); // why "true" is not working
	
	// =====
	// RIP
	// =====
	// test code for map_vertices_to_circle
	// U = V.leftCols(2);
	// for (int ii=0; ii< bdIdx.rows(); ii++){
	// 	U.row(bdIdx(ii)) = VCir.row(ii);
	// }
}

