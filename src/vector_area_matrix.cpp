#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>

typedef Eigen::Triplet<double> tri;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  /*
    Construct the symmetric matrix A which is the vector-area of the mesh given by F
    This is formed through the minimization of angle distortions in the parameterization:

    min ½ ∫ ‖ ∇u - ∇v⟂ ‖² dA
    
    which can be expanded into 

    min ½ ∫ ‖∇u‖² + ‖∇v‖² - ∇u·∇v⟂ dA
    
    min ½ DirichletNRG(u, v) - ∫∇u·∇v⟂ dA

    focusing on the second term, we can see that
    ∇u·∇v⟂ = det( / ∂u/∂x  ∂u/∂y \
                  \ ∂v/∂x  ∂v/∂y / ) 
    which should equal to 1 in order to minimize angle distortions.

    min ½ DirichletNRG(u, v) - ∫1 dA

    Using Stokes, we can get the same result around the ∂S (the boundary of the mesh)

    min ½ DirichletNRG(u, v) - ∮ u·n ds

    min ½ DirichletNRG(u, v) - ∑∫ ½(uᵢ + t(uᵢ - uⱼ))·(uⱼ - uᵢ)⟂/‖uⱼ - uᵢ‖ dt

    min ½ DirichletNRG(u, v) - ∑ ½(uⱼ - uᵢ)·(uⱼ - uᵢ)

    we want an A such that: UᵀAU = ½ ∑ (uⱼ - uᵢ)·(uⱼ - uᵢ)⟂

    where the uᵢ = (uᵢx uᵢy) = U(i, i+n)

    the A matrix needs to be symmetric for our solver, so set A̅ = ½(A + Aᵀ)
   */

  // get the largest boundary loop (the ∂S)
  std::vector<std::vector<int32_t>> boundaries;
  igl::boundary_loop(F, boundaries);

  std::vector<tri> triplets;
  triplets.reserve(4 * boundaries.size());

  // get the number of vertices:
  int32_t n = F.maxCoeff() + 1;

  for(int32_t j = 0; j < boundaries.size(); j++)
  {
    std::vector<int32_t>boundary = boundaries[j];
    for(int32_t i = 0; i < boundary.size(); i++)
    {
      int32_t iOne = boundary[i];
      int32_t iTwo = boundary[(i + 1) % boundary.size()];

      // Construct the A matrix: (and scale by one half)
      triplets.push_back(tri(iOne, n + iTwo,  0.25));
      triplets.push_back(tri(n + iOne, iTwo, -0.25));

      // Construct the Aᵀ matrix: (and scale by one half)
      triplets.push_back(tri(n + iTwo, iOne,  0.25));
      triplets.push_back(tri(iTwo, n + iOne, -0.25));    
    }
  }

  // fully construct the A matrix:
  A.resize(2*n, 2*n);
  A.setFromTriplets(triplets.begin(), triplets.end());
}

