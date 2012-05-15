/* Simple Objects that provide an easier interface for SWIG wrapper generation */

#ifndef SCRIPT_H_
#define SCRIPT_H_
#include <vector>
#include <fstream>
#include "simpledb.h"
class Descriptors;

class Database {
private:
	std::fstream f;
	SimpleDB* db;
	char mode;
	int check_mode(const char* expected=NULL, const char* error=NULL);
public:
	Database();
	virtual ~Database();
	int open(const char* filename, char mode);
	void rewind();
	void close();
	Descriptors next();
	Descriptors get(unsigned int index);
	int store(Descriptors& d);
};

class Descriptors {
private:
	std::vector<unsigned int> descs;
public:
	Descriptors();
	int parse_sdf(const char* sdf);
	int parse_sdfile(const char* sdfile);
	int parse_smiles(const char* smile);
	unsigned int get_descriptor(unsigned int i);
	unsigned int get_len();
	virtual ~Descriptors();
	friend double similarity(Descriptors*, Descriptors*);
	friend class Database;
};

double similarity(Descriptors* d1, Descriptors* d2);
int batch_parse(const char* sdfile, const char* dbfile);


#endif /*SCRIPT_H_*/
