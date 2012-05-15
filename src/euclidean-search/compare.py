#!/usr/bin/env python
# a script for comparing two results: one from knn, one from db_search
import sys

accuracies = []
try:
	f1, f2 = sys.argv[1:3]
	f1 = file(f1)
	f2 = file(f2)
except:
	sys.stderr.write("Usage: %s result.1 result.2\n" % sys.argv[0])
	sys.exit(1)
lines1 = [set([int(ii) for ii in i.split(',')][1:]) for i in f1]
lines2 = [set([int(ii) for ii in i.split(',')][1:]) for i in f2]
for i, s1 in enumerate(lines1):
	s2 = lines2[i]
	accuracies.append( len(s1.intersection(s2)) * 1.0 / len(s2) )

print ",".join([str(i) for i in accuracies])
print sum(accuracies) / len(accuracies)
print min(accuracies) 
ones = 0
ninety_nine = 0
ninety = 0
for i in accuracies:
	if i == 1:
		ones += 1
	if i >= 0.99:
		ninety_nine += 1
	if i > 0.9:
		ninety += 1
	
print ones
print ninety_nine
print ninety
