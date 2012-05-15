#include "desc.h"
#include "simpledb.h"
#include "search.h"
#include "profiling.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>

void usage()
{
	std::cerr << "Usage: db_search [-rhd] database.cdb query.sdf K" << std::endl;
	std::cerr << "Usage: db_search -i[d] database.cdb queryid.list K" << std::endl;
	std::cerr << "Usage: db_search [-s start(included) -e end(included) -d] database.cdb K" << std::endl;
	std::cerr << "       All indices are 1-based" << std::endl;
}

void output(std::vector<std::pair<double, unsigned int> >& results,
		bool details, 
		std::vector<unsigned int> desc_query,
		IndexedDB & db)
{
		if (details) {
			std::vector<unsigned int> desc_db;
			for (unsigned int i = results.size(); i > 0; i --) {
				db.at(results[i-1].second - 1, desc_db);
				std::cout << results[i-1].second << ":"
					<< 1 - similarity(desc_query, desc_db) << " ";
			}
			std::cout << std::endl;
		} else {
			std::cout << (results.end() - 1)->second;
			for (unsigned int i = 0; i < results.size(); i ++)
				std::cout << "," << results[i].second;
			std::cout << std::endl;
		}

}

int main(int argc, char* argv[])
{
	int c;
	int start = -1;
	int end = -1;
	extern int optind;
	bool reverse = false;
	bool details = false;
	bool use_list = false;
	bool use_range = false;
	unsigned int k;
	const char *db_file, *query_file; 
	while ((c = getopt(argc, argv, "irhds:e:")) != -1) {
		switch(c) {
			case 's': start = atoi(optarg); break;
			case 'e': end = atoi(optarg); break;
			case 'i': use_list = true; break;
			case 'r': reverse = true; break;
			case 'd': details = true; break;
			case 'h': usage(); return 1; break;
			case '?': 
								std::cerr << "Unknown options. use -h to see usage" << std::endl;
								return 1;
		}
	}

	argv = argv + optind - 1;

	db_file = argv[1];
	query_file = NULL;

	if (argc - optind == 2 and start >=0 and end > 0 and end > start) {
		/* valid and range run */
		use_range = true;
		k = atoi(argv[2]);
	} else if (argc - optind != 3) {
		/* invalid arguments */
		usage();
		return 1;
	} else {
		/* valid and not range run */
		k = atoi(argv[3]);
		query_file = argv[2];
	}
	
	if (k == 0) {
		std::cerr << "Invalid K supplied. Must be a positive integer. You supplied ";
		if (use_range) std::cout << argv[2];
		else std::cout << argv[3];
		std::cout << std::endl;
		return 1;		
	}
	
	DEBUG_MSG("read database file.")
	Timer t;
	IndexedDB db;
	t.start();
	if (not db.open(db_file)) return 1;
	std::cerr << "Database loaded in " << t.pause() << " seconds" << std::endl;
	t.reset();

	std::fstream ifs;
	if (not use_range) {
		ifs.open(query_file, std::ios::in);
		if (not ifs.good()) {
			std::cerr << "Cannot open query file for reading. File is "
				<< query_file << std::endl;
			return 1;		
		}
	}

	std::vector<unsigned int> desc_query;
	std::vector<std::pair<double, unsigned int> > results;
	if (not use_list and not use_range) {
		std::string sdf; int line = 0;
		int sdf_id = 0;
		while (sdf_iter(ifs, sdf, line)) {
			std::cerr << "processing query " << ++sdf_id << std::endl;
			Molecule* query = new_mol_from_sdf(sdf.c_str());
			if (query == NULL) {
				std::cerr << "SDF is invalid, near line " << line << std::endl;
				continue;
			}
			desc_query.clear();
			calc_desc(*query, desc_query);
			delete query;
			results.clear();
			t.start();
			knn(db, desc_query, k, results, 1, reverse);
			t.pause();

			output(results, details, desc_query, db);
		}
	} else if (not use_range) {
		unsigned int qid;
		while (true) {
			ifs >> qid;
			if (not ifs.good()) break;
			desc_query.clear();
			db.at(qid - 1, desc_query);
			t.start();
			results.clear();
			knn(db, desc_query, k, results, 1, reverse);
			t.pause();

			output(results, details, desc_query, db);
		}
	} else {
		for (int qid = start; qid <= end; qid ++) {
			desc_query.clear();
			db.at(qid - 1, desc_query);
			t.start();
			results.clear();
			knn(db, desc_query, k, results, 1, reverse);
			t.pause();

			output(results, details, desc_query, db);
		}
	}

	std::cerr << "Total search time is " << t.read() << " seconds." << std::endl;
	return 0;
}

// vim:tabstop=2:shiftwidth=2
