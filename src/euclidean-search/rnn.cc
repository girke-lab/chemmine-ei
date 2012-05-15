#include "eucsearch.h"
#include "search.h"
#include <iostream>
#include <algorithm>
#include "profiling.h"

#ifdef _REFINE_
#include "desc.h"
#include "simpledb.h"
#endif

int main(int argc, const char* argv[])
{
	if (argc != 4 and argc != 7) {
		std::cerr << "Usage: " << argv[0]
			<<" db.record query.record r [db.cdb query.sdf r2]" 
			<< std::endl;
		return 1;
	}

	double r = atof(argv[3]);

	Timer t;
	t.start();

#ifdef _REFINE_
	double r2 = 0; if (argc == 7) r2 = atof(argv[6]);
	PreloadedDB db;
	std::fstream ifs;
	if (argc == 7) {
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

	std::cerr << "Data loaded in " << t.pause() << " seconds" << std::endl;
	t.reset();
	t.start();
	Timer t_wo_r;

	for (int q_id = 0; q_id < query_data.size(); q_id ++) {
		t_wo_r.start();
		std::cerr << "processing query " << q_id << std::endl;
		std::vector<unsigned int> hits;
		radius_search(db_coord, query_coord.begin() + q_id * query_dim,
			query_dim, r, hits);
		std::cerr << query_data[q_id] << " hits " << hits.size()
			<< " results before refining." << std::endl;
		std::cout << query_data[q_id];
		/*
		for (int i = 0; i < hits.size(); i ++)
			std::cout << "," << hits[i];
		std::cout << std::endl;
		*/

		t_wo_r.pause();
#ifdef _REFINE_
		if (argc == 7) {
			if (sdf_iter(ifs, sdf, line_cntr) == 0) {
				std::cerr << "Error in refining: no more SDF" << std::endl;
				return 0;
			}
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
			for (int i = 0; i < hits.size(); i ++) {
				db.at(hits[i] - 1, dbcmp);
				double s = similarity(desc_query, dbcmp);
				if (s > 1 - r2)
					std::cout << "," << hits[i];
			}
			std::cout << std::endl;
		}
#endif
	}

	t.pause();
	std::cerr << "Time: " << t_wo_r.read()
		<< " seconds without refinement" << std::endl;
	std::cerr << "Time: " << t.read() << " seconds total in search" << std::endl;

	return 0;
}
// vim:tabstop=2:shiftwidth=2:smartindent
