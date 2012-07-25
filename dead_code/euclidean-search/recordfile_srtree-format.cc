#include "eucsearch.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
/* a record file is structured after format defined in SRTree impl:
 * line 1 : dimensions
 * line 2 : v1:v2:v3:...:vn:(data)
 * line 3: same as line 2
 */

static int _check_f(std::fstream& f) {
	if (not f.good() and not f.eof()) {
		std::cerr << "** bad format" << std::endl;
		return 0;
	}
	return 1;
}

static void _clean_up(int& dim, std::vector<double> coord, std::vector<datatype> data)
{
	dim = 0;
	data.clear();
	coord.clear();
}

static void _read_line(const char* _line, int dim, 
		std::vector<double>& coord, std::vector<datatype>& data)
{
	char* num = NULL;
	char *endp;
	char *line  = new char[strlen(_line) + 1];
	strcpy(line, _line);
	for (int i = 0; i < dim; i ++) {
		if (num == NULL) 
			num = strtok(line, ":");
		else
			num = strtok(NULL, ":");
		if (endp == num) {
			std::cerr << "invalid value: %s; inserting -2 instead"
				<< std::endl;
			coord.push_back(-2);
		}
		coord.push_back(strtod(num, &endp));
	}
	char* last = strtok(NULL, ":");
	int id = atoi(last + 1);
	data.push_back(datatype(id));
}

int read_file(const char* rec_file, int& dim, std::vector<double>& coord, 
		std::vector<datatype>& data, bool quiet)
{
	std::string line;
	dim = 0;
	data.clear();

	std::fstream f(rec_file, std::ios::in);

	if (not f.good()) {
		std::cerr << "** cannot open file for read" << std::endl;
		return 0;
	}

	f >> dim;
	if (_check_f(f) == 0) return 0;
	if (_check_f(f) == 0) return 0;

	int ln_num = 0;
	while (not f.eof()) {
		f >> line;
		if (_check_f(f) == 0) {
			_clean_up(dim, coord, data);
			return 0;
		} else if (not f.eof()) {
			_read_line(line.c_str(), dim, coord, data);
			ln_num ++;
		}
	}	

	if (not quiet)
		std::cout << ln_num << std::endl;

	return 1;
}

#ifdef _TEST_RECORDFILE_CC
int main(int argc, char* argv[])
{
	int dim;
	std::vector<double> coord;
	std::vector<datatype> data;
	read_file(argv[1], dim, coord, data);
	for (int i = 0; i < data.size(); i ++) {
		for (int j = 0; j < dim; j ++) 
			std::cout << coord[i * dim + j] << " ";
		std::cout << ": " << data[i] << std::endl;
	}
		
	return 0;
}
#endif
