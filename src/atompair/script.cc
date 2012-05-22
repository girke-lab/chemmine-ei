#include "script.h"
#include "desc.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "db_build.h"
#include <string.h>


Descriptors::Descriptors()
{
}

Descriptors::~Descriptors()
{
	//std::cout << "descriptors object released" << std::endl;
}

int Descriptors::parse_sdf(const char* sdf)
{
	Molecule* mol = new_mol_from_sdf(sdf);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

int Descriptors::parse_sdfile(const char* sdfile)
{
	Molecule* mol = new_mol_from_sdfile(sdfile);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

int Descriptors::parse_smiles(const char* smiles)
{
	Molecule* mol = new_mol_from_smiles(smiles);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

unsigned int Descriptors::get_descriptor(unsigned int i)
{
	if (i >= descs.size()) return 0;
	return descs[i];
}

unsigned int Descriptors::get_len()
{
	return descs.size();
}

double similarity(Descriptors* d1, Descriptors* d2)
{
	if (d1 == NULL || d2 == NULL) {
		std::cerr << "one or both input compounds are invalid" << std::endl;
		return 0;
	}
	
	return similarity(d1->descs, d2->descs, 1);
}

int batch_parse(const char* sdfile, const char* dbfile)
{
	return batch_sdf_parse(sdfile, dbfile);
}

Database::Database()
{
	mode = 0;
	db = NULL;
}

Database::~Database()
{
	close();
}
void Database::close()
{
	if (f.is_open()) f.close();
	if (db) {
		delete db;
		db = NULL;
	}
}

int Database::open(const char* filename, char _mode)
{
	close();
	if (_mode == 'r' or _mode == 'i' or _mode == 'm') {
		if (_mode == 'r')
			db = new SimpleDB();
		else if (_mode == 'i')
			db = new IndexedDB();
		else
			db = new PreloadedDB();

		if (not db->open(filename)) {
			delete db;
			db = NULL;
			return 0;
		}
	} else if (_mode == 'w') {
		f.open(filename, std::iostream::out);
		if (not(f.is_open() && f.good())) {
			std::cerr << "Cannot open file for writing" << std::endl;
			return 0;
		}
		init_db(f);
	} else if (_mode == 'a') {
		// open a db to check header
		SimpleDB _db;
		if (not _db.open(filename))  return 0;
		_db.close();
		f.open(filename, std::iostream::out | std::iostream::app);
		if (not (f.is_open() && f.good())) {
			std::cerr << "Cannot open file for writing" << std::endl;
			return 0;
		}
	} else {
		std::cerr << "Unknown mode. Only understand r, w, a" << std::endl;
		return 0;
	}
	mode = _mode;
	
	return 1;
}

int Database::check_mode(const char* expected, const char* error)
{
	if (mode == 0) {
		std::cerr << "Database is not opened yet" << std::endl;
		return 0;		
	}
	
	if (expected == NULL) return 1;
	for (unsigned int i = 0; i < strlen(expected); i ++)
		if (expected[i] == mode) return 1;

	if (error == NULL)
		std::cerr << "The requested operation cannot be performed in current mode";
	else
		std::cerr << error;
	std::cerr << std::endl << "allowd modes: " << expected;
	return 0;
}

Descriptors Database::next()
{
	Descriptors d;
	if (not check_mode("rmi", "read not allowed in this mode")) return d;
	
	if (not db->next(d.descs)) d.descs.clear();
	return d;
}

int Database::store(Descriptors& d)
{
	if (not check_mode("wa", "write not allowed in this mode")) return 0;
	
	return serialize(d.descs, f);
}

void Database::rewind()
{
	if (not check_mode("rmi", "rewind not allowed in this mode")) return;
	db->rewind();
}

Descriptors Database::get(unsigned int index)
{
	Descriptors d;
	if (not check_mode("mi", "random access not allowed in this mode")) return d;

	if (mode == 'm') {
		if (not ((PreloadedDB*) db)->at(index, d.descs)) d.descs.clear();
	} else {
		if (not ((IndexedDB*) db)->at(index, d.descs)) d.descs.clear();
	}
	return d;
}
