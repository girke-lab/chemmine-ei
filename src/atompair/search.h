#include "simpledb.h"
#include <queue>
#include <vector>

int knn(SimpleDB& db, const std::vector<unsigned int>& desc, unsigned int k, 
	std::vector<std::pair<double, unsigned int> >& hits,
	unsigned int index_base=1, bool reverse=false);

int limited_knn(PreloadedDB& db, const std::vector<unsigned int>& desc, unsigned
	int k, std::vector<std::pair<double, unsigned int> >& hits,
	std::vector<unsigned int>& candidates, 
	unsigned int index_base=1, bool reverse=false);

int radius_search(SimpleDB& db, const std::vector<unsigned int>& desc, double r, 
	std::vector<unsigned int>& hits,
	unsigned int index_base=1, bool reverse=false);
