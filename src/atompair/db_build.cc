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
		std::cerr << "Cannot open database file for writing. File is " << dbfile << std::endl;
		ifs.close();
		ofs.close();
		return 0;		
	}

	char cids_f[1024] = "";
	strncpy(cids_f, dbfile, 1017);
	if (strlen(dbfile) >= 1017) cids_f[1016] = '\0';
	strcat(cids_f, ".names");
	std::fstream cids_ofs(cids_f, std::fstream::out);
	if (! cids_ofs.good()) {
		std::cerr << "Cannot open database name file for writing. File is "
							<< cids_f << std::endl;
		ifs.close();
		ofs.close();
		cids_ofs.close();
		return 0;		
	}

	ofs.write(HEADER, HEADER_SIZE);
	char intlen = (char) sizeof(int);
	ofs.write(&intlen, 1);
	
	std::string sdf_buf;
	std::vector<unsigned int> descs;
	int line_cntr = 0;
	int cntr = 0; int writer_cntr = 0;
	while (sdf_iter(ifs, sdf_buf, line_cntr)) {
		std::string cid = sdf_buf.substr(0, sdf_buf.find("\n"));
		Molecule * mol = new_mol_from_sdf(sdf_buf.c_str());
		if (mol) {
			DEBUG_MSG("successfully parsed 1 compound");
			calc_desc(*mol, descs);
			DEBUG_VAR(descs.size());
			int io_stat = serialize(descs, ofs);
			if (io_stat == 2) {
				std::cerr << "Empty compound (#" << cntr << " cid:" << cid
									<< ") ignored" << std::endl;
			} else if (io_stat == 0) {
				std::cerr << "Error when writing to database file" << std::endl;
				ifs.close(); ofs.close(); cids_ofs.close();
				return 0;
			} else {
				// write the name to the index file
				cids_ofs << cid << std::endl;
				writer_cntr ++;
			}
			cntr ++;
		} else {
			std::cerr << "Failed when parsing compound ending at line "
				<< line_cntr << std::endl;
		}
		delete mol;
		descs.clear();
	}
	
	if (ifs.eof()) {
		std::cerr << "reaching end of file" << std::endl;
		std::cerr << cntr << " compounds processed, "
							<< writer_cntr << " compounds written." << std::endl;
	} else {
		std::cerr << "Error when reading line " << line_cntr << std::endl;
		ifs.close(); ofs.close(); cids_ofs.close();
		return 0;
	}
	
	ifs.close(); ofs.close(); cids_ofs.close();
	
	// test output file
	std::cerr << "Verifying generated database" << std::endl;
	SimpleDB db;
	if (not db.open(dbfile)) {
		std::cerr << "Error when opening the database file" << std::endl;
		return 0;
	}
	std::vector<unsigned int> desc;
	cntr = 0;
	while (db.next(desc)) cntr ++;
	std::cerr << cntr << " compounds read." << std::endl;
	if (db.status() != SIMPLEDB_END) {
		std::cerr << "Error when verifying the database file" << std::endl;
		db.close();
		return 0;
	}
	db.close();
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
