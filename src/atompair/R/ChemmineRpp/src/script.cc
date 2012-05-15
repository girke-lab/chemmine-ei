#include "script.h"
#include "desc.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "simpledb.h"
#include "db_build.h"

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
}

Database::~Database()
{
	close();
}
void Database::close()
{
	if (f.is_open()) f.close();
}

int Database::open(const char* filename, char _mode)
{
	if (f.is_open()) f.close();
	if (_mode == 'r') {
		f.open(filename, std::iostream::in);
	} else if (_mode == 'w') {
		f.open(filename, std::iostream::out);
	} else if (_mode == 'a') {
		f.open(filename, std::iostream::out | std::iostream::app);
	} else {
		std::cerr << "Unknown mode. Only understand r, w, a" << std::endl;
		return 0;
	}
	mode = _mode;
	
	if (f.is_open() && f.good())
		return 1;
	return 0;
}

Descriptors Database::next()
{
	Descriptors d;
	if (mode == 0) {
		std::cerr << "Database is not opened yet" << std::endl;
		return d;		
	}
	if (mode != 'r') {
		std::cerr << "Can only read database in  'r' mode" << std::endl;
		return d;
	}
	
	int ret = load(d.descs, f);
	if (ret == 0)
		d.descs.clear();
	return d;
}

int Database::store(Descriptors& d)
{
	if (mode == 0) {
		std::cerr << "Database is not opened yet" << std::endl;
		return 0;		
	}
	if (mode != 'w') {
		std::cerr << "Can only read database in  'w' mode" << std::endl;
		return 0;
	}
	
	return serialize(d.descs, f);
}
