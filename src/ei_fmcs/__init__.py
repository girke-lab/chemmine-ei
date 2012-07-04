
import sys
from eutils.sdfiterator import sdf_iter

def write_subset(orig_db,index,subset):
	out = file(subset,"w")
	def writeSdf(sdf):
		out.write(sdf)
	foreach_compound(orig_db,index,writeSdf)
	out.close()


def foreach_compound(db,index,operation):
	f = file(index)

	line = f.readline()
	if not line:
		sys.exit(0)
	next_index = int(line)

	#print "starting"
	for i,current_sdf in enumerate(sdf_iter(db)):
	#	print "i=%d" % i
		if i+1 < next_index:
			continue
		elif i+1 == next_index:
			operation(current_sdf)
			line = f.readline()
			if not line:
				break;
			next_index = int(line)
		else: # somehow we got ahead of the index
			raise StandardError("index does not seem to be sorted. aborting\n")
	#print "done"
	f.close()
