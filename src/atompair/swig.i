%{
#include "script.h"
#include <vector>
%}
class Descriptors{
public:
   Descriptors();
   int parse_sdf(const char* sdf);
   int parse_sdfile(const char* sdfile);
   int parse_smiles(const char* smile);
   unsigned int get_descriptor(unsigned int i);
   unsigned int get_len();
};

double similarity(Descriptors* d1, Descriptors* d2);
int batch_parse(const char* sdfile, const char* dbfile);

class Database {
public:
   Database();
   int open(const char* filename, char mode);
   void close();
   Descriptors next();
   int store(Descriptors& d);
   void rewind();
   Descriptors get(unsigned int index);
};
