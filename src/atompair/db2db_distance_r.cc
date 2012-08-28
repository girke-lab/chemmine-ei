#include "desc.h"
#include "simpledb.h"
#include <R.h>
#include <Rinternals.h>

#include "helpers.h"

extern "C" {
   SEXP db2db_distance_iddb(SEXP dbFile, SEXP iddb1File,SEXP iddb2File);
   SEXP db2db_distance_db(SEXP dbFile, SEXP db2File);
}


SEXP db2db_distance_iddb(SEXP dbFile, SEXP iddb1File,SEXP iddb2File)
{
  IndexedDB db;
  std::vector<unsigned int> iddb1;
  std::vector<unsigned int> iddb2;
  std::vector<unsigned int> cmp1,  cmp2;
  SEXP ans;

  if(openDb(CHAR(STRING_ELT(dbFile,0)),db)) error("could not open db %s\n",dbFile);
  if(openIddb(CHAR(STRING_ELT(iddb1File,0)),iddb1)) error("could not read iddb1 file %s\n",iddb1File);
  if(openIddb(CHAR(STRING_ELT(iddb2File,0)),iddb2)) error("could not read iddb2 file %s\n",iddb2File);

  PROTECT(ans = allocMatrix(REALSXP,iddb1.size(),iddb2.size()));

  for (unsigned int i = 0; i < iddb1.size(); i ++) {
     db.at(iddb1[i], cmp1);
     for (unsigned int j = 0; j < iddb2.size(); j ++) {
        db.at(iddb2[j], cmp2);
        REAL(ans)[i+iddb1.size()*j]= 1 - similarity(cmp1,  cmp2);
     }
  }

  db.close();
  UNPROTECT(1);
  return ans;
}
SEXP db2db_distance_db(SEXP dbFile, SEXP db2File)
{
  IndexedDB db,db2;
  std::vector<unsigned int> cmp1,cmp2;
  SEXP ans;

  if(openDb(CHAR(STRING_ELT(dbFile,0)),db)) error("could not read db file %s\n",dbFile);
  if(openDb(CHAR(STRING_ELT(db2File,0)),db2)) error("could not read db2 file %s\n",db2File);

  PROTECT(ans = allocMatrix(REALSXP,db.size(),db2.size()));
  for (unsigned int i = 0; i < db.size(); i ++) {
     db.at(i, cmp1);
     for (unsigned int j = 0; j < db2.size(); j ++) {
        db2.at(j, cmp2);
        REAL(ans)[i+db.size()*j]= 1 - similarity(cmp1,  cmp2);
     }
  }

  db2.close();
  db.close();

  UNPROTECT(1);
  return ans;
}
