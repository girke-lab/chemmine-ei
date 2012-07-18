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
#include <cassert>

void usage()
{
	std::cerr << "Usage: db_search [-rhd] database.cdb query.sdf K" << std::endl;
	std::cerr << "Usage: db_search -i[d] database.cdb db.list queryid.list K" << std::endl;
	std::cerr << "Usage: db_search [-s start(included) -e end(included) -d] database.cdb K" << std::endl;
	std::cerr << "Usage: db_search [-d] -s start -e end -f database.cdb K [-o index_offset] raw.results"
		<< std::endl;
	std::cerr << "       -r : reverse (find the farthest)" << std::endl;
	std::cerr << "       -i : input is a file containing indices" << std::endl;
	std::cerr << "       -d : show distances" << std::endl;
	std::cerr << "       -f : refinement mode" << std::endl;
	std::cerr << "       -o : offset to add to indices in raw.results. positive."
		<< std::endl;
	std::cerr << "       -c : continue from file. This skips the first few entries in db by looking at a result file previously generated. Note that candidate file is not skipped. Perfect for progressively-generated-and-refined candidate file." << std::endl;
	std::cerr << "       All indices are 1-based" << std::endl;
}

void output(std::vector<std::pair<double, unsigned int> >& results,
		bool details, 
		std::vector<unsigned int> desc_query,
		PreloadedDB & db)
{
		if (details) {
			std::vector<unsigned int> desc_db;
			for (unsigned int i = results.size(); i > 0; i --) {
				db.at(results[i-1].second - 1, desc_db);
				std::cout << results[i-1].second
				<< ":" << 1 - similarity(desc_query, desc_db) << " ";
			}
			std::cout << std::endl;
		} else {
			std::cout << (results.end()-1)->second;
			for (unsigned int i = 0; i < results.size(); i ++)
				std::cout << "," << results[i].second;
			std::cout << std::endl;
		}

}

const unsigned int BUFSIZE = 10485760;

int main(int argc, char* argv[])
{
	int c;
	int start = -1;
	int end = -1;
	int offset = 0;
	extern int optind;
	bool reverse = false;
	bool details = false;
	bool use_list = false;
	bool refinement = false;
	bool use_range = false;
	bool batch = false;
	unsigned int k;
	const char *db_file, *id_file, *query_file, *candidate_file, *batch_file; 
	batch_file = NULL;
	char *linebuf = NULL;
	std::ifstream *resume_f = NULL;
	while ((c = getopt(argc, argv, "b:c:o:firhds:e:")) != -1) {
		switch(c) {
			case 'b': batch = true; batch_file = optarg; break;
			case 'c': resume_f = new std::ifstream(optarg); break;
			case 's': start = atoi(optarg); break;
			case 'e': end = atoi(optarg); break;
			case 'i': use_list = true; break;
			case 'r': reverse = true; break;
			case 'd': details = true; break;
			case 'f': refinement = true; break;
			case 'h': usage(); return 1; break;
			case 'o': offset = atoi(optarg); break;
			case '?': 
								std::cerr << "Unknown options. use -h to see usage" << std::endl;
								return 1;
		}
	}

	if (resume_f) assert(resume_f->good());

	argv = argv + optind - 1;

	db_file = argv[1];
	id_file = NULL;
	query_file = NULL;
	candidate_file = NULL;

	if (not refinement and argc - optind == 2 and start >=0 and end > 0 and end >= start) {
		/* valid and range run */
		use_range = true;
		k = atoi(argv[2]);
	} else if (refinement) {
		if (argc - optind != 3 or start == 0 or end == 0 or end < start) {
			usage();
			return 1;
		}
		if (offset < 0) {
			std::cerr << "wrong offset value. It must be a positive value."
				<< std::endl;
			return 1;
		}
		candidate_file = argv[3];
		use_range = true;
		k = atoi(argv[2]);
		linebuf = new char[BUFSIZE];
	} else if (batch and argc - optind == 2) {
		k = atoi(argv[2]);
	} else if (argc - optind == 4 and use_list) {
		id_file = argv[2];
		k = atoi(argv[4]);
		query_file = argv[3];
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
	PreloadedDB db;
	t.start();
	if (use_list and id_file != NULL) {
		if (not db.open(db_file, id_file)) return 1;
	} else {
		if (not db.open(db_file)) return 1;
	}
	std::cerr << "Database loaded in " << t.pause() << " seconds" << std::endl;
	t.reset();

	PreloadedDB batch_db;
	if (batch) {
		DEBUG_MSG("read batch file.")
		Timer t;
		t.start();
		if (not batch_db.open(batch_file)) return 1;
		std::cerr << "Batch database loaded in " << t.pause() << " seconds" << std::endl;
		t.reset();
	}

	std::fstream ifs;
	if (not batch and not use_range) {
		ifs.open(query_file, std::ios::in);
		if (not ifs.good()) {
			std::cerr << "Cannot open query file for reading. File is "
				<< query_file << std::endl;
			return 1;		
		}
	}

	std::fstream candiate_f;
	if (refinement) {
		candiate_f.open(candidate_file, std::ios::in);
		if (not candiate_f.good()) {
			std::cerr << "Cannot open result file for reading. File is "
				<< candidate_file << std::endl;
			return 1;		
		}
	}

	std::vector<unsigned int> desc_query;
	std::vector<std::pair<double, unsigned int> > results;
	if (batch) {
		desc_query.clear();
		while (batch_db.next(desc_query)) {
			t.start();
			results.clear();
			knn(db, desc_query, k, results, 1, reverse);
			t.pause();

			output(results, details, desc_query, db);
			desc_query.clear();
		}
	} else if (not use_list and not use_range) {
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
			//std::cout << desc_query.size() << std::endl;
			t.start();
			results.clear();
			if (id_file) db.enable_mask();
			knn(db, desc_query, k, results, 1, reverse);
			if (id_file) db.disable_mask();
			t.pause();

			output(results, details, desc_query, db);
		}
	} else {
		unsigned int skipped = 0;
		for (int qid = start; qid <= end; qid ++) {
			desc_query.clear();
			db.at(qid - 1, desc_query);
			t.start();
			results.clear();
			if (refinement) {
				if (resume_f) {
					resume_f->getline(linebuf, BUFSIZE);
					if (resume_f->good()) {skipped ++; t.pause(); continue;}
					else {
						std::cerr << skipped << " lines skipped." << std::endl;
						resume_f->close();
						delete resume_f;
						resume_f = NULL;
					}
				}
				candiate_f.getline(linebuf, BUFSIZE);
				if (candiate_f.eof()) {
					std::cerr << "Reaching EOF when reading result file." << std::endl;
					return 1;
				}
				if (candiate_f.fail()) {
					std::cerr << "I/O error when reading result file. Line too long?"
						<< std::endl;
					return 1;
				}

				// parsing
				unsigned int id;
				char* str;
				std::vector<unsigned int> candidates;
				str = strtok(linebuf, ":");
				do {
					id = atoi(str);
					assert(id);
					candidates.push_back(id + offset);
					str = strtok(NULL, " ");
					assert(str);
					str = strtok(NULL, ":");
				} while (str);
				limited_knn(db, desc_query, k, results, candidates, 1, reverse);
			} else {
				knn(db, desc_query, k, results, 1, reverse);
			}
			t.pause();

			output(results, details, desc_query, db);
		}
	}

	std::cerr << "Total search time is " << t.read() << " seconds." << std::endl;

	if (linebuf) {
		delete linebuf;
		linebuf = NULL;
	}

	return 0;
}

// vim:tabstop=2:shiftwidth=2
