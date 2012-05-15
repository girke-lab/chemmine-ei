#include "eucsearch.h"
#include "search.h"
#include <iostream>
#include <fstream>

// read index pairs from a file, and output similarity
int main(int argc, const char* argv[])
{

	if (argc != 3) {
		std::cerr << "Usage: " << argv[0]
			<<" db.record index-pair.file" 
			<< std::endl;
		return 1;
	}

	int db_dim;
	std::vector<double> db_coord;
	std::vector<datatype> db_data;
	if (not read_file(argv[1], db_dim, db_coord, db_data)) {
		std::cerr << "cannot process file " << argv[1] << std::endl;
	}

	std::fstream ifs(argv[2], std::ios::in);
	if (not ifs.good()) {
		std::cerr << "Cannot read file " << argv[2] << std::endl;
		return 1;
	}
	unsigned int left, right;
	while (ifs.good()) {
		ifs >> left; ifs >> right;
		if (not ifs.good()) break;
		std::cerr << "processing query " << left << ":" << right
			<< "        \r";
		double d = distf(db_coord.begin() + left * db_dim,
				db_coord.begin() + right * db_dim, db_dim);
		std::cout << 1 - d << std::endl;
	}

	return 0;
}

// vim:tabstop=2:shiftwidth=2:smartindent
