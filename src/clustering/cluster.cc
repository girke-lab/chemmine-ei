/*
 * Patrick-Jarvis Clustering
 * Input: a neighbor file, each line is of this format:
 *        entry_id,neighbor_1,neighbor_2,... neighbor_p
 * Output: to STDOUT
 *        cluster_id : entry_id1, entry_id 2, ....
 *        cluster_id : entry_id, ...
 */
#include "debug.h"
#include "DisjointSets.h"
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cassert>
#include <algorithm>

#define LINE_BUF_SIZE 10240

std::map<std::string, int> name_to_id;
std::vector<std::string> names;
std::vector<std::vector<int> > nbr_list;

/* read in a neighbor file */
void static prepare_neighbors(const char* nbr_file, int skip, int p)
{
	// open file
	std::fstream f(nbr_file, std::ios::in);
	if (not f.good()) {
		std::cerr << "**Error in opening file " << nbr_file << std::endl;
		exit(1);
	}

	// read line by line
	char line_buf[LINE_BUF_SIZE];
	for (int i = 0; i < skip; i ++) {
		f.getline(line_buf, LINE_BUF_SIZE);
	}
	int line_cntr = 0;
	names.clear();
	name_to_id.clear();
	while (f.good()) {
		f.getline(line_buf, LINE_BUF_SIZE);
		// parse the first before comma
		char* raw = strtok(line_buf, ",");
		if (raw == NULL) continue;
		DEBUG_VAR(line_cntr);
		DEBUG_VAR(raw);
		names.push_back(std::string(raw));
		name_to_id[std::string(raw)] = line_cntr;
		line_cntr ++;
	}
	
	assert(f.eof());
	f.clear();
	f.seekg(0, std::ios::beg);
	for (int i = 0; i < skip; i ++) {
		f.getline(line_buf, LINE_BUF_SIZE);
	}

	nbr_list.clear();
	line_cntr = 0;
	while (f.good()) {
		f.getline(line_buf, LINE_BUF_SIZE);
		if (f.eof()) break;
		// parse the first before comma
		char* raw = strtok(line_buf, ",");
		if (raw == NULL) break;
		// now parse the rest
		std::vector<int> nbrs;
		while (true) {
			char* raw = strtok(NULL, ",");
			if (raw == NULL) break;
			DEBUG_VAR(raw);
			int nbr = name_to_id[std::string(raw)];
			if (nbr != line_cntr) nbrs.push_back(nbr);
			DEBUG_VAR(nbr);
		}
		// apply p
		int to_erase = nbrs.size() - p;
		if (to_erase > 0)
			nbrs.erase(nbrs.begin(), nbrs.begin() + to_erase);
		// sort nbrs
		std::sort(nbrs.begin(), nbrs.end());
		nbr_list.push_back(nbrs);
		line_cntr ++;
	}
}

void print_neighbors()
{
	for (int i = 0; i < nbr_list.size(); i ++) {
		std::cout << i;
		for (int j = 0; j < nbr_list[i].size(); j ++) {
			std::cout << " " << nbr_list[i][j];
		}
		std::cout << std::endl;
	}
}

/* size of neighbor set intersection. Expect a sorted list */
int nbr_intersect(std::vector<int>& nbrs1, std::vector<int>& nbrs2)
{
	int i = 0, j = 0, intrsct = 0;
	while (i < nbrs1.size() and j < nbrs2.size()) {
		if (nbrs1[i] == nbrs2[j]) {
			intrsct ++; 
			i ++; j ++;
		} else if (nbrs1[i] > nbrs2[j])
			j ++;
		else
			i ++;
	}
	return intrsct;
}

/* initialize the clustering */
DisjointSets s;
void cluster_init() {
	s.AddElements(names.size());
	return;
}

void cluster(int m)
{
	cluster_init();
	for (int i = 0; i < names.size(); i ++) {
		for (int j = 0; j < nbr_list[i].size(); j ++) {
			int nbr = nbr_list[i][j];
			if (s.FindSet(i) == s.FindSet(nbr)) continue;
			// check condition 2
			if (nbr_intersect(nbr_list[i], nbr_list[nbr]) < m)
				continue;
			// merging clusters
			s.Union(s.FindSet(i), s.FindSet(nbr));
		}
	}
}

void print_clusters()
{
	for (int i = 0; i < names.size(); i ++) 
		std::cout << s.FindSet(i) << std::endl;	
}

void print_cluster_stat(int print_pair = 0)
{
	std::map<int, std::vector<int> > lg_cls;
	int *st = new int[names.size()];
	bzero(st, sizeof(int) * names.size());
	assert(st[0] == 0 and st[1] == 0);
	for (int i = 0; i < names.size(); i ++) {
		int set = s.FindSet(i);
		if (st[set] == 0) {
			st[set] = i + 1;
		} else if (st[set] > 0) {
			std::vector<int> _t;
			_t.push_back(st[set] - 1);
			_t.push_back(i);
			st[set] = -1;
			lg_cls[set] = _t;
		} else {
			assert(st[set] == -1);
			lg_cls[set].push_back(i);
		}
	}

	if (print_pair == 0) {
		for (std::map<int, std::vector<int> >::iterator i = lg_cls.begin();
				i != lg_cls.end();
				i ++) {
			std::cout << "cluster " << (*i).first << ":";
			for (int j = 0; j < (*i).second.size(); j ++)
				std::cout << " " << names[(*i).second[j]];
			std::cout << std::endl;
		}
	} else {
		for (std::map<int, std::vector<int> >::iterator i = lg_cls.begin();
				i != lg_cls.end();
				i ++) {
			for (int j = 0; j < (*i).second.size(); j ++) {
				for (int k = j + 1; k < (*i).second.size(); k ++) {
					long _i = (*i).second[j] * names.size() + (*i).second[k];
					long _j = (*i).second[k] * names.size() + (*i).second[j];
					if (_i > _j) 
						std::cout << " " << _j;
					else
						std::cout << " " << _i;
				}
			}
		}
		std::cout << std::endl;
	}
}

/*
 *
 * Modes:
 * 1: print non-singleton clusters and their members
 * 2: for each point, print its cluster
 * other: print pairs, encoded in a single long int
 */
int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Usage: " << argv[0] << " input.nbr skip_line p m [mode]" << std::endl;
		return 1;
	}
	int mode = 0;
	if (argc == 6)
		mode = atoi(argv[5]);
	prepare_neighbors(argv[1], atoi(argv[2]), atoi(argv[3]));
	//print_neighbors();
	cluster(atoi(argv[4]));
	if (mode == 2)
		print_clusters();
	else if (mode == 1)
		print_cluster_stat();
	else
		print_cluster_stat(1);
	return 0;
}
