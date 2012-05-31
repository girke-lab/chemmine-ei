#include "bdbdb.h"
#include <iostream>
#include <vector>
#include <fstream>
#define BUFSIZE 81

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cerr << "Usage: bdb_import in.db in.names out.bdb" << std::endl;
		return 1;
	}

	SimpleDB db;
	if (db.open(argv[1]) == 0) {
		std::cerr << "Cannot open database " << argv[1] << " for reading."
			<< std::endl;
		return 1;
	}

	char linebuf[BUFSIZE];
	std::ifstream ifs(argv[2], std::ios::in);
	if (not ifs.good()) {
		std::cerr << "Cannot open file " << argv[2] << std::endl;
		db.close();
		return 1;
	}

	BDBDB bdb;
	std::vector<unsigned int> cmp;
	if (bdb.open(argv[3], 'c') == 0) {
		std::cerr << "Cannot open database " << argv[3] << " for writing."
			<< std::endl;
		db.close();
		ifs.close();
		return 1;
	}

	std::cout << "importing..." << std::endl;

	while(db.next(cmp)) {
		ifs.getline(linebuf, BUFSIZE);
		if (not ifs.good()) {
			std::cerr << "Cannot open file " << argv[2] << std::endl;
			db.close();
			bdb.close();
			ifs.close();
			return 1;
		}
		std::cout << "writing " << linebuf << "           \r";
		bdb.write(linebuf, cmp);
	}

	db.close();
	bdb.close();
	ifs.close();

	std::cout << "verifying...                    " << std::endl;
	if (bdb.open(argv[3]) == 0) {
		std::cerr << "Cannot open database " << argv[3] << " for reading."
			<< std::endl;
		return 1;
	}
	unsigned int counter = 0;
	std::string buf;
	while(bdb.next(cmp, buf)) {
		counter ++;
		std::cout << buf << "          \r";
	}
	std::cout << counter << " compounds imported" << std::endl;
	return 0;
}
