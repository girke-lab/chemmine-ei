#include "bdbdb.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "desc.h"

int main(int argc, char* argv[])
{
	unsigned char mode;
	bool pass = true;
	if (argc == 3)
		mode = 'a';
	else if (argc == 4)
		mode = 'c';
	else 
		pass = false;

	int n = 0;
	
	if (pass and mode == 'a') {
		n = atoi(argv[2]);
		if (n == 0) 
			pass = false;
	} 

	if (not pass) {
		std::cerr << "Usage: bdb_verify name.bdb count" << std::endl;
		std::cerr << "Usage: bdb_verify name.bdb cid1 cid2" << std::endl;
		return 1;
	}

	BDBDB bdb;

	if (bdb.open(argv[1]) == 0) {
		std::cerr << "Cannot open database " << argv[1] << " for reading."
			<< std::endl;
		return 1;
	}
	
	if (mode == 'a') {
		unsigned int counter = 0;
		std::string buf;
		std::vector<unsigned int> cmp;
		while(bdb.next(cmp, buf) and n-- > 0) {
			counter ++;
			if (counter % 1000 == 0)
				std::cout << counter << "                   \r";
		}
		std::cout << std::endl;
		std::cout << counter << " compounds verified" << std::endl;
	} else if (mode == 'c') {
		std::vector<unsigned int> cmp1, cmp2;
		bdb.read(argv[2], cmp1);
		bdb.read(argv[3], cmp2);
		std::cout << "similarity: " << similarity(cmp1, cmp2) << std::endl;
	}
	return 0;
}
