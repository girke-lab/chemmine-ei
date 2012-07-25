#include "eucsearch.h"
#include "search.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include "profiling.h"

#ifdef _REFINE_
#include "desc.h"
#include "simpledb.h"
#endif

int main(int argc, const char* argv[])
{

	if (argc != 4 and argc != 7) {
		std::cerr << "Usage: " << argv[0]
			<<" db.record query.record k [db.cdb query.sdf k2]" 
			<< std::endl;
		return 1;
	}
	int k = atoi(argv[3]);

	Timer t;
	t.start();

#ifdef _REFINE_
	int k2 = 0;
	PreloadedDB db;
	std::fstream ifs;
	if (argc == 7) {
		k2 = atoi(argv[6]);
		if (not db.open(argv[4])) return 1;

		ifs.open(argv[5], std::fstream::in);
		if (! ifs.good()) {
			std::cerr << "Cannot open SDF file for reading. File is "
				<< argv[5] << std::endl;
			ifs.close();
			return 1;		
		}
	}
	std::string sdf;
	int line_cntr = 0;
#endif

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

	double elapsed = t.pause();
	std::cerr << "Data loaded in " << elapsed << " seconds" << std::endl;
	t.reset();
	t.start();
	Timer t_wo_r;

	for (unsigned int q_id = 0; q_id < query_data.size(); q_id ++) {
		t_wo_r.start();
		std::cerr << "processing query " << q_id << std::endl;
		std::vector<unsigned int> hits;
		knn(db_coord, query_coord.begin() + q_id * query_dim,
			query_dim, k, hits);
		std::cout << query_data[q_id];
		if (argc != 7) {
			for (unsigned int i = 0; i < hits.size(); i ++)
				std::cout << "," << hits[i];
			std::cout << std::endl;
		}
		t_wo_r.pause();

#ifdef _REFINE_
		if (argc == 7) {
			if (sdf_iter(ifs, sdf, line_cntr) == 0) {
				std::cerr << "Error in refining: no more SDF" << std::endl;
				return 0;
			}
			//std::cerr << "refining search ..." << std::endl;
			Molecule* query = new_mol_from_sdf(sdf.c_str());
			if (query == NULL) {
			std::cerr << "SDF is invalid, near line " << line_cntr
				<< ". Will exit."
				<< std::endl;
			return 1;
			}

			std::vector<unsigned int> desc_query;
			calc_desc(*query, desc_query);
			delete query;
			std::vector<unsigned int> dbcmp;
			std::priority_queue<std::pair<double, unsigned int> > q_dist;
			for (int i = 0; i < hits.size(); i ++) {
				db.at(hits[i] - 1, dbcmp);
				double s = similarity(desc_query, dbcmp);
				if (q_dist.size() == k2) {
					if (1-s > q_dist.top().first) continue;
					q_dist.pop();
					q_dist.push(std::pair<double, unsigned int>(1-s, hits[i]));
				} else {
					q_dist.push(std::pair<double, unsigned int>(1-s, hits[i]));
				}
	
			}
			while (not q_dist.empty()) {
				std::cout << "," << q_dist.top().second;
				q_dist.pop();
			}
			std::cout << std::endl;
		}
#endif
	}

	std::cerr << "Time: " << t_wo_r.read() << " seconds without refinement"
		<< std::endl;
	std::cerr << "Time: " << t.pause() << " seconds total in search" << std::endl;
	std::cerr << "Output Format: query,neighbor_n,neighbor_n-1,...neighbor_1" << std::endl;
	return 0;
}

// vim:tabstop=2:shiftwidth=2:smartindent
