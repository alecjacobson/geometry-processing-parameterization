#include "active_set.h"
#include <set>

template <>
const ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Free>() const {
    return m_active_pair;
}
template <>
const ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Dirichlet>() const {
    return m_dirichlet_pair;
}

template <>
ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Free>() {
    return m_active_pair;
}
template <>
ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Dirichlet>() {
    return m_dirichlet_pair;
}


void ActiveSet::updateActive() {
    std::set<Index> constrained;

    auto&& diri = getIndices<ActiveSet::VariableType::Dirichlet>();
    std::copy(diri.begin(),diri.end(),std::inserter(constrained,constrained.end()));
    int prev = 0;



    auto& active_indices = m_active_pair.indices;
    active_indices.clear();
    active_indices.reserve(active_size());

    for(auto&& d: constrained) {
        if(prev < d) 
            for(; prev < d; ++prev) {
                active_indices.push_back(prev);
            }
        ++prev;
    }
    if(prev != m_total_size) {
        for(; prev < m_total_size; ++prev) {
            active_indices.push_back(prev);
        }
    }


    m_active_pair.updateInverse(m_total_size); 


}
void ActiveSet::IndexPair::updateInverse(size_t total_size) {
    inverse = std::vector<Index>(total_size,-1);
    for(int i = 0; i < indices.size(); ++i) {
        inverse[indices[i]] = i;
    }
}



/*
int main(int argc, char* argv[]) {
    std::vector<int> bad;
    int full_size = 20;
    std::vector<int> dir = {0,3,4,10};
  ActiveSet coords(20);

  coords.setVariableType<ActiveSet::VariableType::Dirichlet>(dir);
  coords.updateActive();
  for(auto&& i: coords.getIndices()) {
      std::cout <<i << " ";
  }
  std::cout << std::endl;


  Eigen::SparseMatrix<double> M(full_size,full_size);
    std::vector<Eigen::Triplet<double>> trips;
    for(int i = 0; i < full_size; ++i) {
        trips.emplace_back(i,i,i+1);
    }
    for(int i = 0; i < full_size-1; ++i) {
        trips.emplace_back(i,i+1,-i-1);
        trips.emplace_back(i+1,i,-i-1);
    }

    M.setFromTriplets(trips.begin(),trips.end());
    std::cout << M << std::endl;

    std::cout << "Reduced: " << std::endl;

    std::cout << coords.reduce<ActiveSet::VariableType::Dirichlet>(M) << std::endl;
    std::cout << coords.reduce<ActiveSet::VariableType::Free>(M) << std::endl;
    std::cout << coords.reduceMatrix<ActiveSet::VariableType::Dirichlet,ActiveSet::VariableType::Free>(M) << std::endl;
    std::cout << coords.reduceMatrix<ActiveSet::VariableType::Free,ActiveSet::VariableType::Dirichlet>(M) << std::endl;


}
*/
