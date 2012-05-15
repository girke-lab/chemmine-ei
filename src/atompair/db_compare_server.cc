// based on compare.cc, keeping only db-db-distance mode
#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
#include <signal.h>
#define BUFSIZE 10240

const char* other_db = "/tmp/__db_compare.in";

void on_sig(int status) {return;}
int main(int argc, char *argv[])
{
	char in[BUFSIZE] = "";
	char out[BUFSIZE] = "";
	
	signal(SIGINT, on_sig);
	PreloadedDB db1, db2;
	if (db1.open(argv[1]) == 0) return 1;
	std::cout << "ready" << std::endl;
	std::cout.flush();
	while (true) {
		db1.rewind();
		pause();
		std::ifstream ifs;
		ifs.open(other_db, std::ios::in);
		ifs.getline(in, BUFSIZE);
		if (not ifs.good()) {
			std::cerr << "Input:invalid input" << std::endl;
			continue;
		}
		ifs.getline(out, BUFSIZE);
		if (ifs.fail()) {
			std::cerr << "Input:invalid input" << std::endl;
			continue;
		}
		ifs.close();

		if (strlen(in) >= BUFSIZE - 1 or strlen(out) >= BUFSIZE - 1) {
			std::cerr << "Input:invalid path to input database" << std::endl;
			continue;
		}

		std::ofstream ofs;
		ofs.open(out);
		if (not ofs.good()) {
			std::cerr << "Input:invalid output file specified" << std::endl;
			continue;
		}

		if (db2.open(in) == 0) {
			std::cerr << "Input: error in opening" << std::endl;
			continue;
		}
		std::vector<unsigned int> cmp1;
		std::vector<unsigned int> cmp2;
		db2.at(0, cmp2);
		std::cout << "OK:" << std::endl;
		while (db1.next(cmp1)) {
			ofs << 1 - similarity(cmp1,  cmp2) << " ";
		}
		ofs << std::endl;
		db2.close();
		ofs.close();
	}
	db1.close();
  return 0;
}
