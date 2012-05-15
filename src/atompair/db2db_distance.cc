#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
int main(int argc, char *argv[])
{
  if (argc != 4 and argc != 3) {
    std::cerr << "Usage: " << argv[0] << " chem.db 1.iddb 2.iddb" << std::endl;
    std::cerr << "       " << argv[0] << " chem.db chem2.db" << std::endl;
    std::cerr << "        1.iddb and 2.iddb are <iddb>s. Are indices inside are considered to be 1-based" << std::endl;
    exit(1);
  }

  IndexedDB db;
  if (db.open(argv[1]) == 0) {
    std::cerr << "Cannot load database " << argv[1] << std::endl;
    return 1;
  }
	if (argc == 3) {
		IndexedDB db2;
		if (db2.open(argv[2]) == 0) {
			std::cerr << "Cannot load database " << argv[2] << std::endl;
			return 1;
		}

		std::vector<unsigned int> cmp1;
		std::vector<unsigned int> cmp2;
		for (unsigned int i = 0; i < db.size(); i ++) {
			std::cerr << "processing query " << i + 1 << "\r";
			db.at(i, cmp1);
			for (unsigned int j = 0; j < db2.size(); j ++) {
				db2.at(j, cmp2);
				std::cout << 1 - similarity(cmp1,  cmp2) << " ";
			}
			std::cout << std::endl;
		}

		db2.close();
	} else {
		std::vector<unsigned int> iddb1, iddb2;
		std::fstream f;
		f.open(argv[2], std::iostream::in);
		if (! f.good()) {
			std::cerr << "Cannot load iddb " << argv[2] << std::endl;
			return 1;
		}
		unsigned int i;
		f >> i;
		while (f.good()) {
			if (i == 0) {
				std::cerr << "indices must be 1-based!" << std::endl;
				return 1;
			}
			iddb1.push_back(i - 1);
			f >> i;
		}
		f.close();

		f.open(argv[3], std::iostream::in);
		if (! f.good()) {
			std::cerr << "Cannot load iddb " << argv[3] << std::endl;
			return 1;
		}
		f >> i;
		while (f.good()) {
			if (i == 0) {
				std::cerr << "indices must be 1-based!" << std::endl;
				return 1;
			}
			iddb2.push_back(i - 1);
			f >> i;
		}
		f.close();

		std::vector<unsigned int> cmp1;
		std::vector<unsigned int> cmp2;
		for (unsigned int i = 0; i < iddb1.size(); i ++) {
			std::cerr << "processing query " << i + 1 << "\r";
			db.at(iddb1[i], cmp1);
			for (unsigned int j = 0; j < iddb2.size(); j ++) {
				db.at(iddb2[j], cmp2);
				std::cout << 1 - similarity(cmp1,  cmp2) << " ";
			}
			std::cout << std::endl;
		}
	}
  db.close();

	std::cerr << std::endl;
  return 0;
}
