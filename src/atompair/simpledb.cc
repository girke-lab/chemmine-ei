#include "simpledb.h"
#include "debug.h"
#include <cassert>
#include <fstream>
#include <queue>
#include <string.h>
#include <algorithm>

static int load(std::vector<unsigned int> & descs, std::fstream & ifs);
static int load(unsigned int ** p_buf, std::fstream & ifs);

int init_db(std::fstream & ofs)
{
	// write header
	ofs.write(HEADER, HEADER_SIZE);
	char intlen = (char) sizeof(int);
	ofs.write(&intlen, 1);
	if (ofs.good())
		return 1;
	else {
		std::cerr << "Error when writing database header" << std::endl;
		return 0;
	}
}

int serialize(std::vector<unsigned int> & descs, std::fstream & ofs)
{
	int size = descs.size();
	if (size == 0) {
		DEBUG_MSG("Compound with 0 descriptors. Ignored it in the database")
		return 2;
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
	descs.clear();
	unsigned int * buf = NULL;
	int size;
	while ((size = load(&buf, ifs)) < 0) continue;	//skip deleted compounds
	if (size == 0) return 0;
	copy(buf, buf+size, inserter(descs, descs.begin()));
	
	delete[] buf;
	return 1;
}

int load(unsigned int** p_buf, std::fstream & ifs)
{
	int _size, size;
	ifs.read((char*) &_size, sizeof(int));
	if (ifs.fail()) return 0;
	if (_size < 0) size = - _size;
	else size = _size;
	unsigned int * buf = new unsigned int[size];
	ifs.read((char*) buf, sizeof(unsigned int) * size);
	if (ifs.fail()) return 0;
	
	DEBUG_VAR(size);
	*p_buf = buf;
	return _size;
}

SimpleDB::SimpleDB()
{
	begin = 0;
	ifs = NULL;
}

SimpleDB::~SimpleDB()
{
	close();
}

int SimpleDB::open(const char* path)
{
	if (ifs) {
		std::cerr << "The database object is already opened when you try to open again!" << std::endl;
		return 0;
	}
	ifs = new std::fstream(path, std::fstream::in);
	if (! ifs->good()) {
		std::cerr << "Cannot open database file for reading. File is " << path << std::endl;
		ifs->close();
		return 0;		
	}
	char header[HEADER_SIZE];
	ifs->read(header, HEADER_SIZE);
	if (strcmp(header, HEADER) != 0) {
		std::cerr << "The database file is not valid: wrong header" << std::endl;
		ifs->close();
		return 0;				
	}
	char intsize;
	ifs->read(&intsize, 1);
	if (intsize != (char) sizeof(int)) {
		std::cerr << "Database architecture mismatch: this database is made on a different architecture." << std::endl;
		ifs->close();
		return 0;						
	}

	begin = ifs->tellg();
	return 1;
}

void SimpleDB::close()
{
	if (ifs) ifs->close();
	ifs = NULL;
}

int SimpleDB::next(std::vector<unsigned int> & desc)
{
	if (ifs == NULL) {
		std::cerr << "Reading database when it is not opened!" << std::endl;
		return 0;
	}
	
	int ret = load(desc, *ifs);
	if (ret == 0 && not ifs->eof())
		std::cerr << "Error when reading the database file" << std::endl;
	return ret;
}

void SimpleDB::rewind()
{
	set_offset(begin);
}

std::streampos SimpleDB::get_offset()
{
	if (ifs == NULL) {
		std::cerr << "Getting offset of database when it is not opened!"
				<< std::endl;
		return -1;
	}
	else
		return ifs->tellg();
}

void SimpleDB::set_offset(std::streampos offset)
{
	if (ifs == NULL) {
		std::cerr << "Setting offset of database when it is not opened!"
				<< std::endl;
	}
	else {
		ifs->clear();
		ifs->seekg(offset, std::ios::beg);
	}
}

unsigned int SimpleDB::status()
{
	if (ifs == NULL) return SIMPLEDB_CLOSED;
	if (ifs->eof()) return SIMPLEDB_END;
	if (ifs->good()) return SIMPLEDB_GOOD;
	return SIMPLEDB_BAD;
}

unsigned int SimpleDB::abs_index(unsigned int index)
{
	return index;
}

PreloadedDB::PreloadedDB():use_mask(false)
{
	opened = false;
	index = 0;
}

PreloadedDB::~PreloadedDB()
{
}

int PreloadedDB::open(const char* path, const char* mask)
{
	if (mask) {
		std::fstream *mask_fs = new std::fstream(mask, std::fstream::in);
		if (! mask_fs->good()) {
			std::cerr << "Cannot open ID file for reading. File is " << mask
				<< std::endl;
			mask_fs->close();
			return 0;		
		}
		unsigned int id;
		while (mask_fs->good()) {
			*mask_fs >> id;
			if (mask_fs->good())
				selected.push_back(id - 1);
		}
		std::cerr << "Read " << selected.size() << " lines of IDs. Treat them as"\
			" 1-based" << std::endl;
		mask_fs->close();
		delete mask_fs;
		std::sort(selected.begin(), selected.end());
	}
	if (not SimpleDB::open(path)) return 0;
	unsigned int *buf = NULL;
	int size;
	int cntr = 0;
	DEBUG_MSG("preloading .... ")
	while (true) { 
		size = load(&buf, *ifs);
		if (size != 0) {
			std::pair<int, unsigned int*> _p;
			_p.first = size; _p.second = buf;
			db_content.push_back(_p);
			cntr ++;
		}
		if (not ifs->good()) break;
	}

	SimpleDB::close();
	opened = true;
	index = 0;
	return cntr;
}

void PreloadedDB::close()
{
	if (opened) {
		for (unsigned int i = 0; i < db_content.size(); i ++) 
			if (db_content[i].second) delete db_content[i].second;
		db_content.clear();
		opened = false;
		index = 0;
		selected.clear();
	}
}

void PreloadedDB::enable_mask()
{
	if (selected.size())
		use_mask = true;
	else
		std::cerr << "Enable mask when no compound is selected. Ignored!" <<
		std::endl;
}

void PreloadedDB::disable_mask()
{
	use_mask = false;
}

int PreloadedDB::next(std::vector<unsigned int> & descs)
{
	if (not opened) {
		std::cerr << "Reading database when it is not opened!" << std::endl;
		return 0;
	}
	unsigned int _index = index;
	if (use_mask) {
		if (_index >= selected.size()) {
			DEBUG_MSG("Reaching the end.")
			return 0;
		}
		_index = selected[_index];
	}
	if (_index >= db_content.size()) {
		DEBUG_MSG("Reaching the end.")
		return 0;
	}
	int size = db_content[_index].first;
	unsigned int* buf = db_content[_index].second;
	descs.clear();
	copy(buf, buf+size, inserter(descs, descs.begin()));
	index ++;

	return size;
}

void PreloadedDB::rewind()
{
	index = 0;
}

std::streampos PreloadedDB::get_offset()
{
	return (std::streampos) index;
}

void PreloadedDB::set_offset(std::streampos offset)
{
	unsigned int i = (unsigned int) offset;
	index = i;
}

unsigned int PreloadedDB::status()
{
	if (not opened) return SIMPLEDB_CLOSED;
	if (use_mask)
		if (index == selected.size()) return SIMPLEDB_END;
	else
		if (index == db_content.size()) return SIMPLEDB_END;
	return SIMPLEDB_GOOD;
}

int PreloadedDB::at(unsigned int i, std::vector<unsigned int> & desc)
{
	int backup = index;
	index = i;
	int ret = next(desc);
	index = backup;
	return ret;
}

unsigned int PreloadedDB::size()
{
	if (status() == SIMPLEDB_CLOSED) {
		std::cerr << "Database is closed when you try to retrieve the size"
			<< std::endl;
		return 0;
	}

	if (use_mask)
		return selected.size();
	return db_content.size();
}

unsigned int PreloadedDB::abs_index(unsigned int index)
{
	if (use_mask)
		return selected[index];
	return index;
}

IndexedDB::IndexedDB()
{
}

IndexedDB::~IndexedDB(){}

int IndexedDB::open(const char* path)
{
	db_index.clear();
	if (not SimpleDB::open(path)) return 0;
	/* build the index */
	std::vector<unsigned int> dbcmp;
	do {
	db_index.push_back(get_offset());
	} while (next(dbcmp)); 
	db_index.pop_back(); 	// the last entry will always be the end of file
	rewind();
	return 1;
}

void IndexedDB::close()
{
	db_index.clear();
	return SimpleDB::close();
}

int IndexedDB::at(unsigned int i, std::vector<unsigned int> & desc)
{
	if (i >= db_index.size()) return 0;
	std::streampos temp = get_offset();
	std::streampos p = db_index[i];
	set_offset(p);
	int ret = next(desc);
	set_offset(temp);
	return ret;
}

unsigned int IndexedDB::size()
{
	if (status() == SIMPLEDB_CLOSED) {
		std::cerr << "Database is closed when you try to retrieve the size"
			<< std::endl;
		return 0;
	}

	return db_index.size();
}

// vim:shiftwidth=2:tabstop=2:smartindent 
