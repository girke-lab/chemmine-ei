#!/usr/bin/env python
# generate random pairs
import sys

try:
	limit = int(sys.argv[1])
	lower = int(sys.argv[2])
	upper = int(sys.argv[3])
except:
	sys.stderr.write("Usage: %s count range-lower range-upper\n" % sys.argv[0])
	sys.stderr.write("       Both range ends are included\n")
	sys.exit(1)

import random
r = random.Random()
for i in range(limit):
	left = r.randint(lower, upper)
	right = r.randint(lower, upper)
	sys.stdout.write("%d %d\n" % (left, right))

