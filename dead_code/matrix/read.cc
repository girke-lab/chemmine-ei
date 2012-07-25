#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "profiling.h"

void usage(const char* prg)
{
	std::cerr << prg << " input" << std::endl;
}

void verify(const char* fp)
{
	std::fstream ofs;
	int dim;
	ofs.open(fp, std::ios::in);
	unsigned int h;
	float v;
	for (unsigned int i = 0; i < 3; i ++) {
		ofs.read((char*) &h, sizeof(h));
		std::cout << h << " ";
	}
	dim = h;
	std::cout << std::endl;
	unsigned int i = 0;
	while (ofs.good()) {
		ofs.read((char*) &v, sizeof(v));
		if (not ofs.good()) break;
		std::cout << v << " ";
		i ++;
		if (i % dim == 0)
			std::cout << std::endl;
	}
	ofs.close();
}

int main(int argc, char* argv[])
{
	if (argc != 2) {usage(argv[0]); return 1;}
	std::fstream ifs, ofs;
	ifs.open(argv[1], std::ios::in);
	assert(ifs.good());

	verify(argv[1]);

	return 0;

}
