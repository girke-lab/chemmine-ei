#include <queue>
#include <vector>

double distf(std::vector<double>::iterator dbi,
	std::vector<double>::iterator qi, unsigned int dim, 
	int sqr=0);

int knn(std::vector<double>& db, std::vector<double>::iterator qi,
	unsigned int dim, unsigned int k, std::vector<unsigned int>& hits,
	unsigned int index_base=1);

int radius_search(std::vector<double>& db, std::vector<double>::iterator qi, 
	unsigned int dim, double r, std::vector<unsigned int>& hits,
	unsigned int index_base=1);
