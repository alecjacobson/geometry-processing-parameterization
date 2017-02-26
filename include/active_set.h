#ifndef ACTIVE_SET_H
#define ACTIVE_SET_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class ActiveSet {
    public:
        enum VariableType { Free, Dirichlet };
        using Index = int;
        using IndexList = std::vector<int>;

        struct IndexPair {
            IndexList indices;
            IndexList inverse;
            size_t size() const { return indices.size(); }
            void updateInverse(size_t total_size);
        };


        size_t active_size() const { return m_active_pair.size(); }
        size_t reduced_size() const { return active_size(); }
        template <VariableType VarType>
            size_t var_size() const { return getIndexPair<VarType>().size();}

        void updateActive();

        ActiveSet(int k): m_total_size(k) {}
        template <VariableType VarType>
            void setVariableType(const IndexList& indices);

        template <VariableType VarType=VariableType::Free>
        const IndexPair& getIndexPair() const;
        template <VariableType VarType=VariableType::Free>
        IndexPair& getIndexPair() ;


        template <VariableType VarType=VariableType::Free>
        const IndexList& getIndices() const { return getIndexPair<VarType>().indices; }
        template <VariableType VarType=VariableType::Free>
        const IndexList& getInverseIndices() const { return getIndexPair<VarType>().inverse; }

        template <VariableType VarType=VariableType::Free>
        Index reduce(Index i) const { return getIndices<VarType>()[i]; }
        template <VariableType VarType=VariableType::Free>
        Index expand(Index i) const { return getInverseIndices<VarType>()[i]; }

        template <VariableType VarType, typename Derived, typename RetType=typename Eigen::internal::remove_all<typename Eigen::internal::eval<Derived>::type>::type>
             RetType reduce(const Eigen::MatrixBase<Derived>& M) const;

        template <VariableType VarType, typename Derived, typename RetType=typename Eigen::internal::remove_all<typename Eigen::internal::eval<Derived>::type>::type>
             RetType expand(const Eigen::MatrixBase<Derived>& M) const;


        template <VariableType VarType, typename T>
            Eigen::SparseMatrix<T> reduce(const Eigen::SparseMatrix<T>& A) const;
        template <ActiveSet::VariableType VarType, ActiveSet::VariableType VarType2, typename T>
            Eigen::SparseMatrix<T> reduceMatrix(const Eigen::SparseMatrix<T>& A) const;
    private:
        Index m_total_size;
        IndexPair m_active_pair;
        IndexPair m_dirichlet_pair;
};

template <>
const ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Free>() const;
template <>
const ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Dirichlet>() const;

template <>
ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Free>() ;
template <>
ActiveSet::IndexPair& ActiveSet::getIndexPair<ActiveSet::VariableType::Dirichlet>() ;




template <ActiveSet::VariableType VarType>
void ActiveSet::setVariableType(const IndexList& indices) {
    auto&& pair = getIndexPair<VarType>();
    pair.indices = indices;
    pair.updateInverse(m_total_size);
}

template <ActiveSet::VariableType VarType, typename Derived, typename RetType>
RetType ActiveSet::expand(const Eigen::MatrixBase<Derived>& M) const {
    RetType R = RetType::Zero(m_total_size,M.cols());

    auto&& pair = getIndexPair<VarType>();
    assert(M.rows() == pair.size() );
    for(int i = 0; i < pair.size(); ++i) {

        R.row(pair.indices[i]) = M.row(i);
    }

    return R;
}

template <ActiveSet::VariableType VarType, typename Derived, typename RetType>
RetType ActiveSet::reduce(const Eigen::MatrixBase<Derived>& M) const {
    RetType R = RetType::Zero(active_size(),M.cols());

    assert(M.rows() == m_total_size);

    auto&& pair = getIndexPair<VarType>();
    for(int i = 0; i < pair.size(); ++i) {

        R.row(i) = M.row(pair.indices[i]);
    }

    return R;
}






template <ActiveSet::VariableType VarType, typename T>
Eigen::SparseMatrix<T> ActiveSet::reduce(const Eigen::SparseMatrix<T>& M) const {
    return reduceMatrix<VarType,VarType,T>(M);
}

template <ActiveSet::VariableType VarType, ActiveSet::VariableType VarType2, typename T>
Eigen::SparseMatrix<T> ActiveSet::reduceMatrix(const Eigen::SparseMatrix<T>& M) const {
    auto&& pair = getIndexPair<VarType>();
    auto&& pair2 = getIndexPair<VarType2>();


    assert(M.rows() == m_total_size);
    assert(M.cols() == m_total_size);

    Eigen::SparseMatrix<T> A(pair.size(), pair2.size());
    std::vector<Eigen::Triplet<T>> trips;
    for(int i = 0; i < pair2.size(); ++i) {

        int k = pair2.indices[i];

        for(typename Eigen::SparseMatrix<T>::InnerIterator it(M,k); it; ++it) {
            if(pair.inverse[it.row()] != -1) {
                int r = pair.inverse[it.row()];
                int c = pair2.inverse[it.col()];
                trips.emplace_back(r,c,it.value());
            }
        }
    }

    A.setFromTriplets(trips.begin(),trips.end());
    return A;
}


#endif//ACTIVE_SET_H
