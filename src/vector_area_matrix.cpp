#include "vector_area_matrix.h"
#include "igl/boundary_loop.h"

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
    typedef Eigen::Triplet<double> Tr;
    int V_size = F.maxCoeff()+1;
    Eigen::SparseMatrix<double> At;
    std::vector<Tr> trips;
    std::vector< std::vector<int> > bl;
    
    A.resize(V_size*2,V_size*2);
    igl::boundary_loop(F, bl);
    for (auto i : bl) {
        for (int j = 0; j < i.size(); j++) {
            int ij[] = {i[j], i[(j+1) % i.size()]};
            trips.push_back(Tr(ij[0], ij[1] + V_size, 0.5));
            trips.push_back(Tr(ij[0] + V_size, ij[1], -0.5));
        }
    }
    A.setFromTriplets(trips.begin(), trips.end());
    At = A.transpose();
    A = (A + At) * 0.5;
}
