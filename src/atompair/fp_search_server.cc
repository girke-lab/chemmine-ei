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
#include <signal.h>

void usage()
{
	std::cerr << "Usage: db_isearch database.cdb K" << std::endl;
	std::cerr << "       All indices are 1-based" << std::endl;
}

void on_sig(int status)
{
	std::cerr << "waked up to work: " << status << std::endl;
	return;
}

void output(std::vector<std::pair<double, unsigned int> >& results,
		bool details, 
		std::vector<unsigned int> desc_query,
		PreloadedDB & db)
{
		std::cout << "OK:";
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

const unsigned int BUFSIZE = 10485760;

int main(int argc, char* argv[])
{
	signal(SIGINT, on_sig); 
	unsigned int k;
	const char *db_file; 
	char *linebuf = new char[BUFSIZE];

	db_file = argv[1];
	k = atoi(argv[2]);

	if (argc != 3) {
		usage();
		return 1;
	}
	
	if (k == 0) {
		std::cerr << "Invalid K. Must be a positive integer." << std::endl
			<<" You supplied " << argv[2] << std::endl;
		return 1;		
	}
	
	DEBUG_MSG("read database file.")
	Timer t;
	PreloadedDB db;
	t.start();
	if (not db.open(db_file)) return 1;
	std::cerr << "Database loaded in " << t.pause() << " seconds" << std::endl;
	t.reset();


	while (true) {
		unsigned int _k = k;
		std::cout << ">>";
		std::cout.flush();
		pause();

		// parsing input
		char* query_dbf = "query.fp_cdb";
		char* candidate_f = "candidates.data";

		// open query database
		PreloadedDB query_db;
		if (not query_db.open(query_dbf)) {
			std::cerr << "Cannot open " << query_dbf << "!" << std::endl;
			continue;
		}

		// open candidate
		bool use_candidate = false;
		std::vector<unsigned int> candidates;
		std::fstream candidate_fo;
		if (candidate_f) {
			candidate_fo.open(candidate_f, std::ios::in);
			if (not candidate_fo.good()) {
				std::cerr << "Cannot open result file for reading. File is "
					<< candidate_f << std::endl;
				continue;		
			}

			candidate_fo.getline(linebuf, BUFSIZE);
			if (candidate_fo.fail()) {
				std::cerr << "I/O error when reading result file. Line too long?"
					<< std::endl;
				continue;		
			}
			candidate_fo.close();

			// parsing
			unsigned int id;
			char* str;
			str = strtok(linebuf, ":");
			do {
				if (strcmp(str, "/K") == 0) {
					str = strtok(NULL, " ");
					unsigned int __k = atoi(str);
					if (__k > 0) _k = __k;
					str = strtok(NULL, ":");
				} else {
					id = atoi(str);
					assert(id);
					candidates.push_back(id);
					str = strtok(NULL, " ");
					if (str == NULL) break;
					str = strtok(NULL, ":");
				}
			} while (str);
		}

		if (candidates.size() > 0) 
			use_candidate = true;

		// load query. Here we load only 1
		std::vector<unsigned int> desc_query;
		std::vector<std::pair<double, unsigned int> > results;
		desc_query.clear();
		query_db.at(0, desc_query);
		query_db.close();

		// search
		t.reset();
		t.start();
		results.clear();
		if (use_candidate) 
			limited_knn(db, desc_query, _k, results, candidates, 1, false);
		else
			knn(db, desc_query, _k, results, 1, false);

		// print result
		std::cerr << "Query takes " << t.pause() << " seconds" << std::endl;
		output(results, true, desc_query, db);
	}	

	if (linebuf) {
		delete[] linebuf;
		linebuf = NULL;
	}

	return 0;
}

// vim:tabstop=2:shiftwidth=2
