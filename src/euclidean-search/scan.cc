#include "eucsearch.h"
#include "search.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>

void usage(const char* argv[])
{
	std::cerr << "Usage: " << argv[0]
		<<" db.record query.record <max number of results>" 
		<< std::endl;
}

int main(int argc, const char* argv[])
{

	if (argc != 4) {
		usage(argv);
		return 1;
	}

	int n_results = atoi(argv[3]);
	if (n_results <= 0) {
		usage(argv);
		return 1;
	}

	int db_dim;
	std::vector<double> db_coord;
	std::vector<datatype> db_data;
	if (not read_file(argv[1], db_dim, db_coord, db_data)) {
		std::cerr << "cannot process file " << argv[1] << std::endl;
	}

	int query_dim;
	std::vector<double> query_coord;
	std::vector<datatype> query_data;
	if (not read_file(argv[2], query_dim, query_coord, query_data)) {
		std::cerr << "cannot process file " << argv[2] << std::endl;
		return 1;
	}

	if (query_dim != db_dim) {
		std::cerr << "dimensions do not match: " 
			<< query_dim << " != " << db_dim << std::endl;
		return 1;
	}

	for (size_t q_id = 0; q_id < query_data.size(); q_id ++) {
		std::cerr << "processing query " << q_id << "\r";
		std::vector<std::pair<double, datatype> > results;
		for (size_t db_id = 0; db_id < db_data.size(); db_id ++) {
			double d = distf(db_coord.begin() + db_id * db_dim,
				query_coord.begin() + q_id * query_dim, db_dim);
			std::pair<double, datatype> p(d, db_data[db_id]);
			results.push_back(p);
		}
		sort(results.begin(), results.end());

		if (n_results > results.size()) n_results = results.size();

		for (unsigned int i = 0; i < n_results; i ++)
			std::cout << results[i].second << ":" << results[i].first << " ";
		std::cout << std::endl;
	}

	return 0;
}

// vim:tabstop=2:shiftwidth=2:smartindent
