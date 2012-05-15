#include "simpledb.h"
#include <assert.h>

int serialize(std::vector<unsigned int> & descs, std::fstream & ofs)
{
	int size = descs.size();
	if (size == 0) {
		std::cerr << "Compound with 0 descriptors. Ignored it in the database" << std::endl;
		return 1;
	}
	
	unsigned int * buf = new unsigned int[size];
	std::copy(descs.begin(), descs.end(), buf);
	DEBUG_VAR(size);
	
	ofs.write((char*) &size, sizeof(int));
	if (ofs.fail()) return 0;
	ofs.write((char*) buf, sizeof(unsigned int) * size);
	if (ofs.fail()) return 0;
	delete[] buf;
	return 1;
}

int load(std::vector<unsigned int> & descs, std::fstream & ifs)
{
	unsigned int * buf = NULL;
	int size = load(&buf, ifs);
	if (size == 0) return 0;
	copy(buf, buf+size, inserter(descs, descs.begin()));
	
	delete[] buf;
	return 1;
}

int load(unsigned int** p_buf, std::fstream & ifs)
{
	int size;
	ifs.read((char*) &size, sizeof(int));
	if (ifs.fail()) return 0;
	unsigned int * buf = new unsigned int[size];
	ifs.read((char*) buf, sizeof(unsigned int) * size);
	if (ifs.fail()) return 0;
	
	DEBUG_VAR(size);
	*p_buf = buf;
	return size;
}
