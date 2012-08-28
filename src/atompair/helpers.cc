#include "helpers.h"
using namespace std;

int openDb(const char* dbFile,IndexedDB &db)
{
  if (db.open(dbFile) == 0) {
    std::cerr << "Cannot load database " << dbFile << std::endl;
    return 1;
  }
  return 0;
}

int openIddb(const char* iddbFile, vector<unsigned int> &iddb)
{
  std::fstream f;
  f.open(iddbFile, std::iostream::in);
  if (! f.good()) {
     std::cerr << "Cannot load iddb " << iddbFile << std::endl;
     return 1;
  }
  unsigned int i;
  f >> i;
  while (f.good()) {
     if (i == 0) {
        std::cerr << "indices must be 1-based!" << std::endl;
        return 1;
     }
     iddb.push_back(i - 1);
     f >> i;
  }
  f.close();

  return 0;
}


