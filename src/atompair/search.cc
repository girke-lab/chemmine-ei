#include <iostream>
#include "desc.h"
#include "search.h"
#include <fstream>

int knn(SimpleDB& db, const std::vector<unsigned int>& desc, unsigned int k, 
	std::vector<std::pair<double, unsigned int> >& hits,
	unsigned int index_base, bool reverse)
{
	if (db.status() == SIMPLEDB_CLOSED) {
		std::cerr << "Database is closed" << std::endl;
		return 0;
	}

	std::streampos mark = db.get_offset();
	db.rewind();

	std::vector<unsigned int> db_cmp;
	std::priority_queue<std::pair<double, unsigned int> > q_dist;
	unsigned int i = 0;
	double s, v;
	while (db.next(db_cmp)) {
		s = similarity(db_cmp, const_cast<std::vector<unsigned int>&>(desc), 1);
		if (reverse) v = s;
		else v = 1 - s;
		/* keep the smallest v's in the queue */
		if (q_dist.size() == k) {
			if (v < q_dist.top().first) {
				q_dist.pop();
				q_dist.push(std::pair<double, unsigned int>(v,
					db.abs_index(i) + index_base));
			}
		} else {
			q_dist.push(std::pair<double, unsigned int>(v,
				db.abs_index(i) + index_base));
		}
		i ++; 
	}
	if (db.status() != SIMPLEDB_END) {
		return 0;
	}
	db.set_offset(mark);

	/* copy result */
	while (! q_dist.empty()) {
		hits.push_back(q_dist.top());
		q_dist.pop();
	}
	return 1;
}

int limited_knn(PreloadedDB& db, const std::vector<unsigned int>& desc,
	unsigned int k, std::vector<std::pair<double, unsigned int> >& hits,
	std::vector<unsigned int> &candidates, unsigned int index_base, bool reverse)
{
	if (db.status() == SIMPLEDB_CLOSED) {
		std::cerr << "Database is closed" << std::endl;
		return 0;
	}

	std::vector<unsigned int> db_cmp;
	std::priority_queue<std::pair<double, unsigned int> > q_dist;
	double s, v;
	for (unsigned int j = 0; j < candidates.size(); j ++) {
		db.at(candidates[j] - index_base, db_cmp);
		s = similarity(db_cmp, const_cast<std::vector<unsigned int>&>(desc), 1);
		if (reverse) v = s;
		else v = 1 - s;
		/* keep the smallest v's in the queue */
		if (q_dist.size() == k) {
			if (v < q_dist.top().first) {
				q_dist.pop();
				q_dist.push(std::pair<double, unsigned int>(v, candidates[j]));
			}
		} else {
			q_dist.push(std::pair<double, unsigned int>(v, candidates[j]));
		}
	}
	std::cerr<<std::endl;

	/* copy result */
	while (! q_dist.empty()) {
		hits.push_back(q_dist.top());
		q_dist.pop();
	}
	return 1;
}

int radius_search(SimpleDB& db, const std::vector<unsigned int>& desc,
	double r, std::vector<unsigned int>& hits,
	unsigned int index_base, bool reverse)
{
	if (db.status() == SIMPLEDB_CLOSED) {
		std::cerr << "Database is closed" << std::endl;
		return 0;
	}

	std::streampos mark = db.get_offset();
	db.rewind();

	std::vector<unsigned int> db_cmp;
	unsigned int i = index_base;
	while (db.next(db_cmp)) {
		double s = similarity(db_cmp, const_cast<std::vector<unsigned int>&>(desc),
				1);
		if (reverse) {
			if (s < 1 - r) hits.push_back(i);
		} else {
			if (s > 1 - r) hits.push_back(i);
		}
		i ++;
	}
	if (db.status() != SIMPLEDB_END) 
		return 0;
	db.set_offset(mark);
	
	return 1;
}
