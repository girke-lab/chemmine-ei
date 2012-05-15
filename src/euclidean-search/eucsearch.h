typedef int datatype;
#include <vector>
int read_file(const char* rec_file, int& dim, std::vector<double>& coord, 
		std::vector<datatype>& data, bool quiet=false);
