#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " chem.db index.file sub.cdb" << std::endl;
    std::cerr << "        indices in index.file are 1-based" << std::endl;
    exit(1);
  }

  IndexedDB db;
  if (db.open(argv[1]) == 0) {
    std::cerr << "Cannot load database " << argv[1] << std::endl;
    return 1;
  }

	std::vector<unsigned int> ids;
	std::fstream f;
	f.open(argv[2], std::iostream::in);
	if (! f.good()) {
		std::cerr << "Cannot load index file " << argv[2] << std::endl;
		return 1;
	}
	unsigned int i;
	f >> i;
	while (f.good()) {
		if (i == 0) {
			std::cerr << "indices must be 1-based!" << std::endl;
			return 1;
		}
		ids.push_back(i - 1);
		f >> i;
	}
	f.close();

	f.open(argv[3], std::iostream::out);
	if (! f.good()) {
		std::cerr << "Cannot open " << argv[3] << " for write." << std::endl;
		return 1;
	} else {
		std::cerr << "Opening " << argv[3] << " for write." << std::endl;
	}

	f.write(HEADER, HEADER_SIZE);
	char intsize = (char) sizeof(int);
	f.write(&intsize, 1);
	std::vector<unsigned int> cmp;
	for (unsigned int i = 0; i < ids.size(); i ++) {
		db.at(ids[i], cmp);
		serialize(cmp, f);
	}
  db.close();

	std::cerr << "Wrote " << ids.size() << " entries." << std::endl;
  return 0;
}
