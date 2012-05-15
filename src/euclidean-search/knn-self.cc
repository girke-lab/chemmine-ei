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

	if (argc != 5 and argc != 7) {
		std::cerr << "Usage: " << argv[0]
			<<" db.record start end k [db.cdb k2]" << std::endl
			<< " indexes are 1-based" << std::endl;
		return 1;
	}
	int s = atoi(argv[2]);
	int e = atoi(argv[3]);
	int k = atoi(argv[4]);

	Timer t;
	t.start();

#ifdef _REFINE_
	int k2 = 0;
	IndexedDB db;
	std::fstream ifs;
	if (argc == 7) {
		k2 = atoi(argv[6]);
		if (not db.open(argv[5])) return 1;

	}
#endif

	int db_dim;
	std::vector<double> db_coord;
	std::vector<datatype> db_data;
	if (not read_file(argv[1], db_dim, db_coord, db_data, true)) {
		std::cerr << "cannot process file " << argv[1] << std::endl;
	}

	double elapsed = t.pause();
	std::cerr << "Data loaded in " << elapsed << " seconds" << std::endl;
	t.reset();
	t.start();
	Timer t_wo_r;

	for (unsigned int q_id = s - 1; q_id < db_data.size() and q_id < e; q_id ++) {
		t_wo_r.start();
		if (q_id % 100 == 0)
			std::cerr << "processing query " << q_id << "          \r";
		std::vector<unsigned int> hits;
		knn(db_coord, db_coord.begin() + q_id * db_dim,
			db_dim, k, hits);
		if (argc != 7) {
			for (unsigned int i = 0; i < hits.size(); i ++)
				std::cout << "," << hits[i];
			std::cout << std::endl;
		}
		t_wo_r.pause();

#ifdef _REFINE_
		if (argc == 7) {
			std::vector<unsigned int> desc_query;
			db.at(q_id, desc_query);
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
			std::vector<std::pair<double, unsigned int> > r(k2);
			int idx = 0;
			while (not q_dist.empty()) {
				r[idx ++] = q_dist.top();
				q_dist.pop();
			}
			for (int i = 0; i < k2; i ++)
				std::cout << r[k2-i-1].second << ":" << r[k2-i-1].first << " ";
			std::cout << std::endl;
		}
#endif
	}

	std::cerr << "Time: " << t_wo_r.read() << " seconds without refinement"
		<< std::endl;
	std::cerr << "Time: " << t.pause() << " seconds total in search" << std::endl;
	return 0;
}

// vim:tabstop=2:shiftwidth=2:smartindent
