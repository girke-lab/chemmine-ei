// based on icoord.cc
#include <stdlib.h>
#include <fstream>
#include <iostream>
#ifdef _TIMER_
#include "profiling.h"
#endif
#include "solver.h"
#include <string.h>
#include <signal.h>
#define LINE_BUF_LIMIT 1048576
using namespace std;
const char* exchange_fp = "/tmp/__coord_server.in";

void on_sig(int status)
{
	return;
}

int read_head(ifstream &ifs, int &k, int &m, double* &p)
{
	// read 1st line
	ifs >> k >> m;
	if (ifs.fail() or ifs.eof()) {
		cerr << "(EE) " << "I/O failed when reading k and m" << endl;
		return 0;
	}
	cerr << "(II) " << "k = " << k << " m = " << m << endl;

	// alloc space
	p = new double[m*k];

	// read p
	for (int i = 0; i < m; i ++)
		for (int j = 0; j < k; j ++) {
			ifs >> p[i*k + j];
			if (ifs.fail() or ifs.eof()) {
				cerr << "(EE) " << "I/O failed when reading p, i = " << i << " j = " << j << endl;
				return 0;
			}
		}

	cerr<< "(II) " << "Finish reading p" << endl;
	return 1;
}


/*
char* f_readline(std::ifstream & ifs)
{
	static char linebuf[LINE_BUF_LIMIT];
	ifs.getline(linebuf, LINE_BUF_LIMIT);
	if (strlen(linebuf) == LINE_BUF_LIMIT - 1) {
		std::cerr << "Input: line way too long" << std::endl;
		return NULL;
	}

	if (ifs.eof()) {
		return NULL;
	}

	if (ifs.fail()) {
		std::cerr << "Input: I/O error" << std::endl;
		return NULL;
	}

	return linebuf;
}


int read_line(int size, double *d)
{
	double _d;
	char* line = NULL;

	std::ifstream ifs;
	ifs.open(exchange_fp, std::ios::in);
	line = f_readline(ifs);
	ifs.close();

	ifs.open(line, std::ios::in);
	line = f_readline(ifs);
	ifs.close();

	if (line == NULL) {
		cerr << "Input: File I/O failed when reading first line" << endl;
		return 0;
	}

	char* str = strtok(line, " ");
	for (int i = 0; i < size; i ++) {
		if (str == NULL) {
			cerr << "Input: File I/O failed when reading d, i = " << i << endl;
			return 0;
		}
		_d = atof(str);
		d[i] = _d;
		if (i < size - 1)
			str = strtok(NULL, " ");
	}

	return 1;
}
*/
int read_line(istream &ifs, int size, double *d)
{
	double _d;
	for (int i = 0; i < size; i ++) {
		ifs >> _d; d[i] = _d;
		
		if (ifs.fail() or ifs.eof()) {
			if (ifs.eof() and i == 0) {
				//cerr << endl << "(II) " << "reaching the end, now stop and exit" << endl;
				return 0;
			}else if(ifs.bad()){
				cerr<<endl<<"(II) stream went bad"<<endl;
			} else {
				cerr << "(EE) " << "I/O failed when reading d, i = " << i << endl;
				throw i;
			}
		}
		
	}
	//cerr<<"done reading"<<endl;

	return 1;
}
int main(int argc, char* argv[])
{
#ifdef _TIMER_
	Timer t;
#endif
	if (argc != 2) {
		cerr << "Usage: " << argv[0] << " puzzle.file" << endl;
		return 1;
	}
	// open file
	ifstream ifs(argv[1], ifstream::in);
	if (not ifs.good()) {
		cerr << "(EE) " << "Cannot open " << argv[1] << endl;
		return 1;
	}
	cerr<< "(II) " << "file opened" << endl;
	int k, m;
	double *p;
	if (read_head(ifs, k, m, p) == 0) {
		cerr << "Error in reading header information" << endl;
		return 1;
	}
	ifs.close();

	Solver s(k, m, M, p);
	delete p;

	int line_id = 0;
	double *d = new double[m];
	doublereal *x = new doublereal[k];
	while (true) {
		if (read_line(cin,m, d) == 0) 
			break;
#ifdef _TIMER_
		t.start();
#endif
		s.optim(x, d);
#ifdef _TIMER_
		t.pause();
#endif
		for (int i = 0; i < k; i ++)
			cout << x[i] << " ";
		cout << endl;
	}
#ifdef _TIMER_
	cerr << "Time: " << t.read() << " seconds" << endl;
#endif

	delete d, x;
	return 0;
}

extern "C" {
	int MAIN__(int argc, char* argv[]){
		main(argc, argv);
	}
}

