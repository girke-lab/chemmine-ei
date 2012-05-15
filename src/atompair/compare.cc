#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
int main(int argc, char *argv[])
{
  if ((argc != 5 && argc != 4) || (
    strcmp(argv[1], "cmp-cmp") != 0 &&
    strcmp(argv[1], "db-cmp") != 0 &&
    strcmp(argv[1], "db-db-distance") != 0 &&
    strcmp(argv[1], "db-db") != 0 &&
    strcmp(argv[1], "bound") != 0 &&
    strcmp(argv[1], "dbself") != 0  
  )) {
    std::cerr << "Usage: " << argv[0] << " cmp-cmp compound1.sdf compound2.sdf" << std::endl;
    std::cerr << "Usage: " << argv[0] << " db-cmp db.cdb compound2.sdf" << std::endl;
    std::cerr << "Usage: " << argv[0] << " db-db db.cdb db2.cdb" << std::endl;
    std::cerr << "Usage: " << argv[0] << " db-db-distance db.cdb db2.cdb" << std::endl;
    std::cerr << "Usage: " << argv[0] << " dbself db.cdb pair.file (indices are 1-based)" << std::endl;
    std::cerr << "Usage: " << argv[0] << " bound db.cdb compound1.sdf compoun2.sdf" << std::endl;
    exit(1);
  }

  if (strcmp(argv[1], "bound") == 0) {
    Molecule* mol1 = new_mol_from_sdfile(argv[3]);
    Molecule* mol2 = new_mol_from_sdfile(argv[4]);
    
    if (mol1 == NULL) {
      std::cerr << "Cannot read SDF file " << argv[3] << std::endl;
      return 1;
    }
    if (mol2 == NULL) {
      std::cerr << "Cannot read SDF file " << argv[4] << std::endl;
      return 1;
    }
    
		std::vector<unsigned int> cmp1;
		std::vector<unsigned int> cmp2;
		calc_desc(*mol1, cmp1);
		calc_desc(*mol2, cmp2);
    std::cout << 1 - similarity(cmp1,  cmp2) << std::endl;
    IndexedDB db;
    if (db.open(argv[2]) == 0) return 1;
    std::vector<unsigned int> dbcmp;
		double lower = 0, upper = 1;
		while (db.next(dbcmp)) {
			double dist1 = 1 - similarity(cmp1, dbcmp);
			double dist2 = 1 - similarity(cmp2, dbcmp);
			if (dist1 - dist2 > lower) lower = dist1 - dist2;
			if (dist2 - dist1 > lower) lower = dist2 - dist1;
			if (dist1 + dist2 < upper) upper = dist1 + dist2;
		}
    std::cout << "[" << lower << "," << upper << "]" << std::endl;
	} else if (strcmp(argv[1], "dbself") == 0) {
		std::fstream ifs(argv[3], std::ios::in);
		if (not ifs.good()) {
			std::cerr << "Cannot read file " << argv[3] << std::endl;
			return 1;
		}
    IndexedDB db;
    if (db.open(argv[2]) == 0) return 1;
		unsigned int left, right;
    std::vector<unsigned int> cmp1;
    std::vector<unsigned int> cmp2;
		while (ifs.good()) {
			ifs >> left; ifs >> right;
			if (not ifs.good()) break;
			std::cerr << "processing query " << left << ":" << right << "     \r";
			db.at(left - 1, cmp1); db.at(right - 1, cmp2);
			std::cout << similarity(cmp1,  cmp2) << " " << std::endl;
		}
	} else if (strcmp(argv[1], "cmp-cmp") == 0) {

    Molecule* mol1 = new_mol_from_sdfile(argv[2]);
    Molecule* mol2 = new_mol_from_sdfile(argv[3]);
    
    if (mol1 == NULL) {
      std::cerr << "Cannot read SDF file " << argv[2] << std::endl;
      return 1;
    }
    if (mol2 == NULL) {
      std::cerr << "Cannot read SDF file " << argv[3] << std::endl;
      return 1;
    }
    
     //std::cout << similarity(mol1,  mol2) << std::endl;
     
     std::vector<unsigned int> v1, v2;
     calc_desc(*mol1, v1);
     calc_desc(*mol2, v2);
     std::cout << similarity(v1,  v2, 1) << std::endl;
    
     delete mol1; delete mol2;
  } else if (strcmp(argv[1], "db-cmp") == 0) {
    SimpleDB db;
    if (db.open(argv[2]) == 0) return 1;
    Molecule* mol2 = new_mol_from_sdfile(argv[3]);
    if (mol2 == NULL) {
      std::cerr << "Cannot read SDF file " << argv[3] << std::endl;
      return 1;
    }
    std::vector<unsigned int> desc;
    calc_desc(*mol2, desc);
    delete mol2;
    std::vector<unsigned int> dbcmp;
    while (db.next(dbcmp)) {
      std::cout << similarity(dbcmp,  desc) << std::endl;
    }
    db.close();
  } else {
		int mode = 0;
		if (strcmp(argv[1], "db-db") == 0) {
			mode = 0;
		} else if (strcmp(argv[1], "db-db-distance") == 0) {
			mode = 1;
		} else {
			std::cerr << "Unknown mode: " << argv[1] << std::endl;
			return 1;
		}
    IndexedDB db1, db2;
    if (db1.open(argv[2]) == 0) return 1;
    if (db2.open(argv[3]) == 0) return 1;
    std::vector<unsigned int> cmp1;
    std::vector<unsigned int> cmp2;
		unsigned cntr = 0;
    while (db1.next(cmp1)) {
			std::cerr << "processing query " << ++cntr << "\r";
      db2.rewind();
      while (db2.next(cmp2)) {
				if (mode == 0)
					std::cout << similarity(cmp1,  cmp2) << " ";
				else
					std::cout << 1 - similarity(cmp1,  cmp2) << " ";
      }
      std::cout << std::endl;
    }
    db1.close();
    db2.close();
		std::cerr << std::endl;
  }
  return 0;
}
