#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
#include "helpers.h"


int db_subset(char* dbFile, char* iddbFile, char* outputFile)
{
   IndexedDB db;
   std::vector<unsigned int> ids;
   if(openDb(dbFile,db)) return 1;
   if(openIddb(iddbFile,ids)) return 1;


   std::fstream f;
	f.open(outputFile, std::iostream::out);
	if (! f.good()) {
		std::cerr << "Cannot open " << outputFile << " for write." << std::endl;
		return -1;
	}

	f.write(HEADER, HEADER_SIZE);
	char intsize = (char) sizeof(int);
	f.write(&intsize, 1);
	std::vector<unsigned int> cmp;
	for (unsigned int i = 0; i < ids.size(); i ++) {
		db.at(ids[i], cmp);
		serialize(cmp, f);
	}
   db.close();

   return ids.size();
}

#ifndef NO_MAIN
int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " chem.db index.file sub.cdb" << std::endl;
    std::cerr << "        indices in index.file are 1-based" << std::endl;
    exit(1);
  }
  int count = db_subset(argv[1],argv[2],argv[3]);
  if(count < 0)
     return -count;
  std::cerr << "Wrote " << count << " entries." << std::endl;
  return 0;
}
#endif
