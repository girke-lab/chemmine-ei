#include "eucsearch.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define BUF_ENTRIES 10000
/* a record file is a binary data file of the following format
 * unsigned int size_of_float; unsigned int num_of_entries; unsigned int dim;
 * float; float; float; ...
 */

int split_file(const char* rec_file, unsigned int n_parts)
{
	unsigned int h[3];

	std::fstream f(rec_file, std::ios::in);

	if (not f.good()) {
		std::cerr << "** cannot open file for read" << std::endl;
		return 0;
	}

	f.read((char*)h, sizeof(h));
	std::cerr << "opening database: float has " << h[0] << " bytes. "
		<< h[1] << " entries. "
		<< h[2] << " dimensions. "
		<< std::endl;
	unsigned int total = h[1];

	unsigned int per_unit = h[1] / n_parts;
	unsigned int offset = 0;
	float *buf = new float[BUF_ENTRIES * h[2]];
	for (unsigned int i = 0; i < n_parts; i ++) {
		unsigned int n_entries = per_unit;
		if (i == n_parts - 1) n_entries = total - per_unit * (n_parts - 1);
		std::string of_name(rec_file);
		of_name += ".";
		char id[8];
		sprintf(id, "%d", i);
		of_name += id;
		std::cerr << "writing " << n_entries << " entries to "
			<< of_name << std::endl;
		std::fstream of(of_name.c_str(), std::ios::out);
		h[1] = (unsigned int) n_entries;
		of.write((char*) h, sizeof(h));
		while (n_entries) {
			std::cerr << ".";
			//std::cerr << "n_entries = " << n_entries << std::endl;
			unsigned int entries_to_read = BUF_ENTRIES;
			if (n_entries < BUF_ENTRIES) entries_to_read = n_entries;
			size_t buf_size = ((size_t) entries_to_read) * h[2] * sizeof(float);
			//std::cerr << "reading " << buf_size << " bytes." << std::endl;
			f.read((char*) buf, buf_size);
			// adjust BUF_ENTRIES when it can read all entries in one read
			assert(f.gcount() == buf_size);
			of.write((char*) buf, buf_size);
			n_entries -= entries_to_read;
		}
		of.close();
		std::cerr << std::endl;
	}
	delete buf;
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "Usage: split matrix.file n_parts" << std::endl;
		return 1;
	}
	unsigned int n = atoi(argv[2]);
	if (n == 0) {
		std::cerr << "Usage: split matrix.file n_parts" << std::endl;
		return 1;
	}

	split_file(argv[1], n);
		
	return 0;
}
