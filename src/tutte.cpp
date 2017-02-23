#include "tutte.h"
#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "edges.h"
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <iostream>

using namespace Eigen;


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	VectorXi b;
	igl::boundary_loop(F, b);

	MatrixXd U_b;
	igl::map_vertices_to_circle(V, b, U_b);

	MatrixXi E = edges(F);

	int n = V.rows(), f = F.rows(), m = E.rows(), nb = b.rows();

	std::vector<bool> isBoundary(n, false);
	std::map<int, int> vToB;
	
	for (int i = 0; i < nb; ++i)
	{
		isBoundary[b(i)] = true;
		vToB[b(i)] = i;
	}

	SparseMatrix<double> A;
	VectorXd B(2*m);
	B.setZero();

	std::vector<Triplet<double>> A_val;

	int row = 0, v1, v2;

	for (int i = 0; i < m; ++i)
	{
		v1 = E(i, 0);
		v2 = E(i, 1);

		if (!isBoundary[v1] && !isBoundary[v2])
		{
			A_val.push_back({ row, v1, 1 });
			A_val.push_back({ row, v2, -1 });
			A_val.push_back({ row + 1, v1 + n, 1 });
			A_val.push_back({ row + 1, v2 + n, -1 });

			row += 2;
		}
		else if (!isBoundary[v1] && isBoundary[v2])
		{
			A_val.push_back({ row, v1, 1 });
			A_val.push_back({ row + 1, v1 + n, 1 });

			B(row) = U_b(vToB[v2], 0);
			B(row + 1) = U_b(vToB[v2], 1);

			row += 2;
		}
		else if (!isBoundary[v2] && isBoundary[v1])
		{
			A_val.push_back({ row, v2, 1 });
			A_val.push_back({ row + 1, v2 + n, 1 });

			B(row) = U_b(vToB[v1], 0);
			B(row + 1) = U_b(vToB[v1], 1);

			row += 2;
		}
	}

	B = B.topRows(row);
	A.resize(row, 2 * n);
	A.setFromTriplets(A_val.begin(), A_val.end());

	std::cout << "Matrix constructed" << std::endl;
	
	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver(A);
	std::cout << "QR decomposition completed" << std::endl;
	
	VectorXd U_in = solver.solve(B);
	std::cout << "Linear system solved" << std::endl;

	U.resize(n, 2);

	for (int i = 0; i < n; ++i)
		if (!isBoundary[i])
		{
			U(i, 0) = U_in(i);
			U(i, 1) = U_in(i + n);
		}

	for (int i = 0; i < nb; ++i)
		U.row(b(i)) = U_b.row(i);
}