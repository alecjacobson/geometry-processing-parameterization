#include "edges.h"
#include <set>

Eigen::MatrixXi edges(const Eigen::MatrixXi &F)
{
  Eigen::MatrixXi E;
  // ADD YOUR CODE HERE

  std::set<std::pair<int, int>> edgeSet;

  for (int i = 0; i<F.rows(); ++i)
  {
	  if (F(i, 0) < F(i, 1))
		  edgeSet.insert(std::make_pair(F(i, 0), F(i, 1)));
	  else
		  edgeSet.insert(std::make_pair(F(i, 1), F(i, 0)));

	  if (F(i, 1) < F(i, 2))
		  edgeSet.insert(std::make_pair(F(i, 1), F(i, 2)));
	  else
		  edgeSet.insert(std::make_pair(F(i, 2), F(i, 1)));

	  if (F(i, 2) < F(i, 0))
		  edgeSet.insert(std::make_pair(F(i, 2), F(i, 0)));
	  else
		  edgeSet.insert(std::make_pair(F(i, 0), F(i, 2)));
  }

  E.resize(edgeSet.size(), 2);

  int count = 0;
  for (std::set<std::pair<int, int> >::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
  {
	  E(count, 0) = it->first;
	  E(count, 1) = it->second;
	  count++;
  }

  return E;
}
