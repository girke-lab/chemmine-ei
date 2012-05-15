#ifndef SIMPLEDB_H_
#define SIMPLEDB_H_
#include "debug.h"
#include <vector>
#include <fstream>
int init_db(std::fstream & ofs);
int serialize(std::vector<unsigned int> & descs, std::fstream & ofs);

#define SIMPLEDB_END  1
#define SIMPLEDB_GOOD  2
#define SIMPLEDB_BAD  0
#define SIMPLEDB_CLOSED 3

class PreloadedDB;
class SimpleDB
{
	private:
		std::streampos begin;
		std::fstream * ifs;
	public:
		SimpleDB();
		virtual ~SimpleDB();
		virtual int open(const char* path);
		virtual void close();
		virtual int next(std::vector<unsigned int> & desc);
		virtual void rewind();
		virtual std::streampos get_offset();
		virtual void set_offset(std::streampos offset);
		virtual unsigned int status();
		virtual unsigned int abs_index(unsigned int);
		friend class PreloadedDB;
};

class PreloadedDB : public SimpleDB
{
	private:
		std::vector<std::pair<int, unsigned int*> > db_content;
		bool opened;
		unsigned int index;
		std::vector<unsigned int> selected;
		bool use_mask;
	public:
		PreloadedDB();
		virtual ~PreloadedDB();
		virtual int open(const char* path, const char* mask=NULL);
		virtual void close();
		virtual int next(std::vector<unsigned int> & desc);
		virtual void rewind();
		virtual std::streampos get_offset();
		virtual void set_offset(std::streampos offset);
		virtual unsigned int status();
		virtual int at(unsigned int i, std::vector<unsigned int> & desc);
		virtual unsigned int size();
		virtual unsigned int abs_index(unsigned int);
		virtual void enable_mask();
		virtual void disable_mask();
};

class IndexedDB : public SimpleDB
{
	private:
		std::vector<std::streampos> db_index;
	public:
		IndexedDB();
		virtual ~IndexedDB();
		virtual int open(const char* path);
		virtual void close();
		virtual int at(unsigned int i, std::vector<unsigned int> & desc);
		virtual unsigned int size();
};

#define HEADER_SIZE 16
#define HEADER "DESCDB"

#endif /*SIMPLEDB_H_*/
