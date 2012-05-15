#include "desc.h"
#include "search.h"
#include "simpledb.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <queue>
#include <stdlib.h>
#include <cassert>

int main(int argc, char* argv[])
{
	if (argc < 3) {
		std::cerr << "Usage: db_search_self database.cdb n [offset limit]" << std::endl;
		return 1;
	}
	unsigned int offset = 0;
	unsigned int limit = 0;
	if (argc > 3) {
		assert(argc == 5);
		offset = atoi(argv[3]);
		limit = atoi(argv[4]);
	}
	
	unsigned int n = atoi(argv[2]);
	if (n == 0) {
		std::cerr << "Invalid n supplied. Must be a positive integer. You supplied " << argv[2] << std::endl;
		return 1;		
	}
	
	PreloadedDB db;
	if (not db.open(argv[1])) return 1;
	unsigned int db_size = db.size();

	std::vector<unsigned int> db_cmp_q;
	std::vector<unsigned int> db_cmp;
	for (unsigned int i = offset; i < db_size; i ++) {
		if (limit and i >= offset + limit) break;
		std::cerr << "processing query " << i+1 << std::endl;
		db.at(i, db_cmp_q);
		std::vector<std::pair<double, unsigned int> > results;
		knn(db, db_cmp_q, n, results);
		std::cout << i + 1;
		for (unsigned int j = 0; j < results.size(); j ++)
			std::cout << "," << results[j].second;
		std::cout << std::endl;
	}
	return 0;
}
