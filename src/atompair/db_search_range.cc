#include "desc.h"
#include "simpledb.h"
#include "search.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <queue>
#include <stdlib.h>

void usage()
{
	std::cerr << "Usage: db_search [-rhd] database.cdb query.sdf radius" << std::endl;
}

int main(int argc, char* argv[])
{
	int c;
	extern int optind;
	bool reverse = false;
	bool details = false;
	while ((c = getopt(argc, argv, "rhd")) != -1) {
		switch(c) {
			case 'r': reverse = true; break;
			case 'd': details = true; break;
			case 'h': usage(); return 1; break;
			case '?': 
								std::cerr << "Unknown options. use -h to see usage" << std::endl;
								return 1;
		}
	}
	if (argc - optind != 3) {
		usage();
		return 1;
	}
	
	argv = argv + optind - 1;
	
	double r = atof(argv[3]);
	if (r == 0) {
		std::cerr << "Invalid n supplied. You supplied " << argv[3] << std::endl;
		return 1;		
	}
	
	IndexedDB db;
	if (not db.open(argv[1])) return 1;
	
	std::string sdf; int line = 0;
	std::fstream ifs(argv[2], std::ios::in);
	if (not ifs.good()) {
		std::cerr << "Cannot open query SDF file for reading. File is "
			<< argv[2] << std::endl;
		return 1;		
	}
	while (sdf_iter(ifs, sdf, line)) {
		Molecule* query = new_mol_from_sdf(sdf.c_str());
		if (query == NULL) {
			std::cerr << "SDF is invalid, near line " << line << std::endl;
			continue;
		}
		std::vector<unsigned int> desc_query;
		calc_desc(*query, desc_query);
		delete query;
		
		std::vector<unsigned int> rst;

		radius_search(db, desc_query, r, rst, 1, reverse);
	
		for (unsigned int i = 0; i < rst.size(); i++) {
			std::cout << rst[i] << " "; 
		}
		std::cout << std::endl;
	}
	return 0;
}
