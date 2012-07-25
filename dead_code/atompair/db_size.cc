#include "desc.h"
#include "simpledb.h"
#include "search.h"
#include "profiling.h"
#include <iostream>
#include <fstream>

void usage()
{
	std::cerr << "Usage: db_search database.cdb" << std::endl;
}

int main(int argc, char* argv[])
{
	if (argc != 2) {
		usage();
		return 1;
	}
	const char* db_file = argv[1];

	Timer t;
	IndexedDB db;
	t.start();
	if (not db.open(db_file)) return 1;
	std::cerr << "Database loaded in " << t.pause() << " seconds" << std::endl;
	std::cout << db.size() << std::endl;
	return 0;
}

// vim:tabstop=2:shiftwidth=2
