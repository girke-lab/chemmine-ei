#include <db_cxx.h>
#include "bdbdb.h"
#include <cassert>
#include "debug.h"
#include "errno.h"
#include <string.h>

BDBDB::BDBDB()
{
	db = NULL;
	cursor = NULL;
	_mode = '_';
}

BDBDB::~BDBDB()
{
	delete db;
	db = NULL;
}

int BDBDB::open(const char* path, const unsigned char _mode)
{
	this->_mode = _mode;
	if (db != NULL) {
		std::cerr << "The database object is already opened when you try to open again!" << std::endl;
		return 0;
	}

	db = new Db(NULL, 0);

	u_int32_t mode;
	if (_mode == 'r')  mode = DB_RDONLY;
	else if (_mode == 'c') mode = DB_CREATE | DB_EXCL;
	else {
		std::cerr << "Unknown mode: " << _mode << std::endl;
		return 0;
	}

	bool pass = false;
	int _errno = 0;
	try{
		db->open(NULL, path, NULL, DB_HASH, mode, 0);
		pass = true;
	} catch(DbException &e) {
		_errno = e.get_errno();
	} catch(std::exception &e) {
	}
	if (not pass) {
		std::cerr << "Cannot open Berkeley DB located at " << path 
			<< ". (requested open mode: " << _mode << ")" << std::endl;
		if (_errno == EEXIST) 
			std::cerr << "Database under the same name exists." << std::endl;
		else if (_errno == ENOENT)
			std::cerr << "Database does not exist." << std::endl;
		else if (_errno == EINVAL)
			std::cerr << "Database is invalid or was not created for hash access."
				<< std::endl;
		return 0;
	}

	char _key[7]; strcpy(_key, "HEADER");
	Dbt key(_key, strlen(_key) + 1);
	Dbt value;
	pass = false;
	if (_mode == 'r') {
		try {
			value.set_flags(DB_DBT_MALLOC);
			assert(db->get(NULL, &key, &value, 0) == 0);
			assert(strcmp(HEADER, (const char*)value.get_data()) == 0);
			pass = true;
		}
		catch(DbException &e) {}
		catch(std::exception &e) {}
	} else if (_mode == 'c') {
		char header[HEADER_SIZE] = HEADER;
		value.set_data(header);
		value.set_size(HEADER_SIZE);
		try {
			assert( db->put(NULL, &key, &value, 0) == 0);
			pass = true;
		}
		catch(DbException &e) {}
		catch(std::exception &e) {}
	}

	if (not pass) {
		std::cerr << "Database initialization failed!" << std::endl;
		try {
			db->close(0);
		}
		catch(DbException &e) {}
		catch(std::exception &e) {}
		return 0;
	}

	pass = false;
	try {
		db->cursor(NULL, &cursor, 0);
		pass = true;
	}
	catch(DbException &e) {}
	catch(std::exception &e) {}

	if (not pass) {
		std::cerr << "Cursor creation failed!" << std::endl;
		try {
			db->close(0);
		}
		catch(DbException &e) {}
		catch(std::exception &e) {}
		return 0;
	}

	return 1;
}

void BDBDB::close()
{
	if (cursor != NULL)
		cursor->close();
	cursor = NULL;

	try {
		db->close(0);
		delete db;
		db = NULL;
	}
	catch(DbException &e) {}
	catch(std::exception &e) {}
}

int BDBDB::next(std::vector<unsigned int> & desc)
{
	std::string c;
	return next(desc, c);
}
int BDBDB::next(std::vector<unsigned int> & desc, std::string& name)
{
	char _name[128];
	int ret = _next(desc, _name);
	if (strcmp(_name, "HEADER") == 0)
		ret = _next(desc, _name);
	name.replace(0, strlen(_name), _name);
	return ret;
}

int BDBDB::_next(std::vector<unsigned int> & desc, char* name)
{
	if (cursor == NULL) {
		std::cerr << "Reading database when it is not opened!" << std::endl;
		return 0;
	}

	Dbt key;
	Dbt value;
	try {
		int ret = cursor->get(&key, &value, DB_NEXT);
		if (ret != 0) {
			assert(ret == DB_NOTFOUND);
			return 0;
		} else {
			if (name) {
				int len = key.get_size();
				strncpy(name, (char*) key.get_data(), len);
			}
			void* buf = value.get_data();
			unsigned int len = value.get_size() / sizeof(int);
			const int * _buf = (const int*) buf;
			std::copy(_buf, _buf + len, std::inserter(desc, desc.begin()));
			return 1;
		}
	} catch(DbException &e) {
		db->err(e.get_errno(), "Error!");
	} catch(std::exception &e) {
		db->errx("Error! %s", e.what());
	}

	std::cerr << "Error when reading the database file" << std::endl;
	return 0;
}

void BDBDB::rewind()
{
	cursor->close();
	db->cursor(NULL, &cursor, 0);
}

std::streampos BDBDB::get_offset()
{
	std::cerr << "get_offset not implemented for Berkeley DB-style DB"
		<< std::endl;
	assert(0);
	return 0;
}

void BDBDB::set_offset(std::streampos offset)
{
	std::cerr << "set_offset not implemented for Berkeley DB-style DB"
		<< std::endl;
	assert(0);
}

unsigned int BDBDB::status()
{
	if (cursor == NULL) return SIMPLEDB_CLOSED;

	Dbc* last;
	cursor->dup(&last, DB_POSITION);
	Dbt key, value;
	int ret = last->get(&key, &value, DB_NEXT);
	if (ret == DB_NOTFOUND) return SIMPLEDB_END;
	last->close();

	try {
		u_int32_t flag;
		assert(db->get_open_flags(&flag) == 0);
		return SIMPLEDB_GOOD;
	}
	catch(DbException &e) {}
	catch(std::exception &e) {}

	return SIMPLEDB_BAD;
}

int BDBDB::write(const char* name, std::vector<unsigned int> & descs)
{
	if (cursor == NULL) {
		std::cerr << "Database is not open" << std::endl;
		return 0;
	}

	if (_mode != 'c') {
		std::cerr << "Cannot write to database in the current open mode."
			<< std::endl;
		return 0;
	}

	Dbt key((void*)name, strlen(name) + 1);
	int size = descs.size();
	if (size == 0) {
		DEBUG_MSG("Compound with 0 descriptors. Ignored it in the database")
		return 2;
	}
	
	unsigned int * buf = new unsigned int[size];
	std::copy(descs.begin(), descs.end(), buf);
	DEBUG_VAR(size);
	Dbt value(buf, sizeof(unsigned int) * size);
	try {
		db->put(0, &key, &value, 0);
		return 1;
	} catch(DbException &e) {
		std::cerr << "Database error when writing." << std::endl;
	} catch(std::exception &e) {
		std::cerr << "General error when writing to database." << std::endl;
	}
	return 0;
}

int BDBDB::read(const char* name, std::vector<unsigned int> & descs)
{
	if (cursor == NULL) {
		std::cerr << "Database is not open" << std::endl;
		return 0;
	}

	if (_mode != 'r') {
		std::cerr << "Cannot read database in the current open mode."
			<< std::endl;
		return 0;
	}

	Dbt key((void*)name, strlen(name) + 1);
	Dbt value;
	try {
		db->get(0, &key, &value, 0);
	} catch(DbException &e) {
		std::cerr << "Database error when writing." << std::endl;
		return 0;
	} catch(std::exception &e) {
		std::cerr << "General error when writing to database." << std::endl;
		return 0;
	}
	void* buf = value.get_data();
	unsigned int len = value.get_size() / sizeof(int);
	const int * _buf = (const int*) buf;
	std::copy(_buf, _buf + len, std::inserter(descs, descs.begin()));
	return 1;
}
