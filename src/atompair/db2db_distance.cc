#include "desc.h"
#include "simpledb.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>
#include "helpers.h"
using namespace std;


int db2db_distance(char *dbFile, char *iddb1File,char *iddb2File,std::ostream &ofs)
{
  IndexedDB db;
  std::vector<unsigned int> iddb1, iddb2;

  if(openDb(dbFile,db)) return 1;
  if(openIddb(iddb1File,iddb1)) return 1;
  if(openIddb(iddb2File,iddb2)) return 1;

  std::vector<unsigned int> cmp1;
  std::vector<unsigned int> cmp2;
  for (unsigned int i = 0; i < iddb1.size(); i ++) {
     std::cerr << "processing query " << i + 1 << "\r";
     db.at(iddb1[i], cmp1);
     for (unsigned int j = 0; j < iddb2.size(); j ++) {
        db.at(iddb2[j], cmp2);
        ofs << 1 - similarity(cmp1,  cmp2) << " ";
     }
     ofs << std::endl;
  }

  db.close();
  std::cerr << std::endl;
  return 0;
}
int db2db_distance(char *dbFile, char *db2File,std::ostream &ofs)
{
  IndexedDB db;
  IndexedDB db2;

  if(openDb(dbFile,db)) return 1;
  if(openDb(db2File,db2)) return 1;

  std::vector<unsigned int> cmp1;
  std::vector<unsigned int> cmp2;
  for (unsigned int i = 0; i < db.size(); i ++) {
     std::cerr << "processing query " << i + 1 << "\r";
     db.at(i, cmp1);
     for (unsigned int j = 0; j < db2.size(); j ++) {
        db2.at(j, cmp2);
        ofs << 1 - similarity(cmp1,  cmp2) << " ";
     }
     ofs << std::endl;
  }

  db2.close();
  db.close();

  std::cerr << std::endl;
  return 0;
}
int db2db_distance2file(char *dbFile, char *iddb1File,char *iddb2File,char *outfile)
{
   std::fstream ofs;
   ofs.open(outfile,std::ios::out);
   int ret=db2db_distance(dbFile,iddb1File,iddb2File,ofs);
   ofs.close();
   return ret;
}
int db2db_distance2file(char *dbFile, char *db2File,char *outfile)
{
   std::fstream ofs;
   ofs.open(outfile,std::ios::out);
   int ret = db2db_distance(dbFile,db2File,ofs);
   ofs.close();
   return ret;
}

#ifndef NO_MAIN
int main(int argc, char *argv[])
{
  if (argc != 4 and argc != 3) {
    std::cerr << "Usage: " << argv[0] << " chem.db 1.iddb 2.iddb" << std::endl;
    std::cerr << "       " << argv[0] << " chem.db chem2.db" << std::endl;
    std::cerr << "        1.iddb and 2.iddb are <iddb>s. All indices inside are considered to be 1-based" << std::endl;
    exit(1);
  }

  if (argc == 4) 
      db2db_distance(argv[1],argv[2],argv[3],std::cout);
  else
      db2db_distance(argv[1],argv[2],std::cout);
  return 0;
}
#endif
