#include "desc.h"
#include "simpledb.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include "db_build.h"

int batch_sdf_parse(const char* sdfile, const char* dbfile)
{
	std::fstream ifs(sdfile, std::fstream::in);
	if (! ifs.good()) {
		std::cerr << "Cannot open SDF file for reading. File is " << sdfile << std::endl;
		ifs.close();
		return 0;		
	}

	std::fstream ofs(dbfile, std::fstream::out);
	if (! ofs.good()) {
		std::cerr << "Cannot open SDF file for writing. File is " << dbfile << std::endl;
		ifs.close();
		ofs.close();
		return 0;		
	}

	char header[HEADER_SIZE] = "DESCDB";
	ofs.write(header, HEADER_SIZE);
	char intlen = (char) sizeof(int);
	ofs.write(&intlen, 1);
	
	std::string sdf_buf;
	std::vector<unsigned int> descs;
	char line[82];
	line[82] = '\0';
	int line_cntr = 0;
	ifs.getline(line, 81);
	line_cntr ++;
	while (ifs.good()) {
		if (strlen(line) == 82) {
			std::cerr << "Line exceeds 80 characters when reading line " << line_cntr << std::endl;
			ifs.close(); ofs.close();
			return 0;
		}
		sdf_buf += line;
		sdf_buf += '\n';
		if (strcmp(line, "$$$$") == 0) {
			Molecule * mol = new_mol_from_sdf(sdf_buf.c_str());
			if (mol) {
				DEBUG_MSG("successfully parsed 1 compound");
				calc_desc(*mol, descs);
				DEBUG_VAR(descs.size());
				serialize(descs, ofs);
			}
			delete mol;
			sdf_buf.clear();
			descs.clear();
		}
		ifs.getline(line, 81);
		line_cntr ++;
	}
	
	if (ifs.eof()) {
		std::cerr << "reaching end of file" << std::endl;
	} else {
		std::cerr << "Error when reading line " << line_cntr << std::endl;
		ifs.close(); ofs.close();
		return 0;
	}
	
	ifs.close(); ofs.close();
	
	// test output file
	std::cerr << "Verifying generated database" << std::endl;
	std::fstream ifs2(dbfile, std::fstream::in);
	std::vector<unsigned int> descs2;
	char header2[HEADER_SIZE];
	ifs2.read(header2, HEADER_SIZE);
	char intlen2;
	ifs2.read(&intlen2, 1);
	if (strcmp(header, header2) != 0) {
		std::cerr << "Error when verifying the database file: headers mismatch." << std::endl;
		return 0;		
	}
	if (intlen2 != intlen) {
		std::cerr << "Error when verifying the database file: integer size mismatch." << std::endl;
		return 0;				
	}
	while (true) {
		descs2.clear();
		if (load(descs2, ifs2) == 0) break;
	}
	if (!ifs2.eof()) {
		std::cerr << "Error when verifying the database file" << std::endl;
		ifs2.close();
		return 0;
	}
	ifs2.close();
	return 1;
}

#ifdef _BUILD_DB_BUILD_
int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "Usage: db_build library.sdf database.cdb" << std::endl;
		return 1;
	}
	
	return 1 - batch_sdf_parse(argv[1], argv[2]);
}
#endif
