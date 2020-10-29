#include "tutte.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void tutte(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		Eigen::MatrixXd & U) {
	// Replace with your code
	// U = V.leftCols(2);

	int vnum = V.rows();
	int fnum = F.rows()

	// calculate L, here, w is 1
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// map the boundary into circle
	// we have some known and some unknown
	Eigen::VectorXd bounds;
	igl::boundary_loop(F, bounds);

	// these points will be maped into circles
	Eigen::MatrixXd UVknown;
	igl::map_vertices_to_circle(V, bounds, UVknown);

	/*
	 *   // MIN_QUAD_WITH_FIXED Minimize a quadratic energy of the form
	 //
	 // trace( 0.5*Z'*A*Z + Z'*B + constant )
	 //
	 // subject to
	 //
	 //   Z(known,:) = Y, and
	 //   Aeq*Z = Beq
	 //   T  should be a eigen matrix primitive type like int or double
	 // Inputs:
	 //   A  n by n matrix of quadratic coefficients
	 //   known list of indices to known rows in Z
	 //   Y  list of fixed values corresponding to known rows in Z
	 //   Aeq  m by n list of linear equality constraint coefficients
	 //   pd flag specifying whether A(unknown,unknown) is positive definite
	 // Outputs:
	 //   data  factorization struct with all necessary information to solve
	 //     using min_quad_with_fixed_solve
	 // Returns true on success, false on error
	 */
	// based on these knows, we solve unknows by minimize trace(UtLU)
	/*
	 *   const Eigen::SparseMatrix<T>& A,
	 const Eigen::MatrixBase<DerivedB> & B,
	 const Eigen::MatrixBase<Derivedknown> & known,
	 const Eigen::MatrixBase<DerivedY> & Y,
	 const Eigen::SparseMatrix<T>& Aeq,
	 const Eigen::MatrixBase<DerivedBeq> & Beq,
	 const bool pd,
	 Eigen::PlainObjectBase<DerivedZ> & Z
	 */

	Eigen::SparseMatrix<double> Aeq;
	igl::min_quad_with_fixed_data<double> data;
	min_quad_with_fixed_precompute(L, bounds, Aeq, false, data);

	Eigen::VectorXd B = Eigen::VectorXd::Zero(vnum);
	Eigen::VectorXd Beq = Eigen::VectorXd::Zero(vnum);
	min_quad_with_fixed_solve(data, B, UVknown, Beq, U);
}

