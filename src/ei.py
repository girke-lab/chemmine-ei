#!/usr/bin/env python
from logging import info, warning, error, debug, critical, root, NOTSET
import sys
import os
import random
from tempfile import NamedTemporaryFile as NTF
import glob
from eutils import OS_Runner, STORED, DISCARDED, getConfig

#sys.path.append('eutils')
#sys.path.append('eutils/fp')
from fpcdb import create_db
from coord import CoordinateSolver
from fpdbcompare import DBComparer
from time import time
from stat import ST_SIZE
from traceback import print_exc
from lshsearch import LSHSearcher
from refineserver import Refiner


root.setLevel(NOTSET)
os_run = OS_Runner()

BINDIR = "" #os.path.join(BASEDIR, 'bin')
DB2DB_DISTANCE = os.path.join(BINDIR, "ei-db2db_distance")
DB_SUBSET = os.path.join(BINDIR,"ei-db_subset")
DB_BUILDER = os.path.join(BINDIR,"ei-db_builder")
EUCSEARCHTOOL = os.path.join(BINDIR, "ei-euclid_search")
EVALUATOR = os.path.join(BINDIR, "ei-evaluator")
COORD_TO_BINARY = os.path.join(BINDIR, "ei-bin_formatter")
INDEXED_SEARCH_EVALUATOR = os.path.join(BINDIR, "ei-comparesearch")
COORDTOOL = os.path.join(BINDIR, "ei-coord")
INDEXED_SEARCH = os.path.join(BINDIR, "ei-isearch")

execfile(getConfig())

BASEDIR =  os.path.abspath(".")
DATADIR = os.path.join(BASEDIR, 'data')
CDB = os.path.join(DATADIR, 'chem.db')
IDDB = os.path.join(DATADIR, 'main.iddb')
TEST_QUERIES = os.path.join(DATADIR, 'test_query.iddb')
CHEMICAL_SEARCH_RESULTS = os.path.join(DATADIR, 'chemical-search.results.gz')
MAX_EUCSEARCH_RESULTS = 50000

INDEXED_SEARCH = INDEXED_SEARCH + " " + lsh_param + " -D %s -C " + CDB + " < " + TEST_QUERIES

class _Cdbsize(object):
	"""a callable that returns the database size and remembers it"""
	def __init__(self):
		self._cdbsize = None
	def __call__(self):
		if self._cdbsize is None:
			info("checking database size...")
			self._cdbsize = \
				int(os_run("wc -l %s" % (IDDB,), stdout=STORED)[1].split()[0])
		return self._cdbsize
cdbsize = _Cdbsize()



def subset_db(n, cdb=IDDB, outfile="ref.iddb", base=1):
	"""build ref cdb by randomly subsetting large cdb"""
	queries = None
	if os.path.exists(TEST_QUERIES):
		queries = set([int(i) for i in file(TEST_QUERIES)])
	if not queries:
		info("No test queries exists. Will create one with %d queries" % 
			n_sample_queries)
		from random import sample
		queries = sample(range(base, cdbsize() + base), n_sample_queries)
		queries.sort()
		f = file(TEST_QUERIES, 'w')
		for q in queries: f.write("%d\n" % q)
		f.close()
		queries = set(queries)
		error("!!!!! PLEASE GENERATE THE CHEMICAL SEARCH RESULT !!!!!")
		error("!!!!! EXAMPLE:                                   !!!!!")
		error("      ei-db_search -id chem.db main.iddb test_query.iddb 50000 | gzip > chemical-search.results.gz")
		error(" ")
		error("Please press any key to continue or CTRL-C to exit")
		sys.stdin.read()
	all = set(range(base, cdbsize() + base))
	assert queries.issubset(all)
	ref = random.sample(all.difference(queries), n)
	ref.sort()

	try:
		f = file(outfile, 'w')
		src = file(cdb)
		for i, line in enumerate(src):
			if not ref: break
			if i + base == ref[0]:
				f.write(line)
				del ref[0]
	except:
		sys.stderr.write("error in subsetting database\n")
		raise

	f.close()
	src.close()
	return outfile

def dist_mat(dbfile, outfile=None):
	"""calculate distance matrix"""
	if not outfile: outfile = dbfile + ".distmat"
	if os.path.exists(outfile):
		debug('reusing small distance matrix file: %s' % outfile)
		return outfile
	cmd = "%s %s %s %s > %s" % (DB2DB_DISTANCE, CDB, dbfile, dbfile,
		outfile)
	os_run(cmd, msg='cannot build distance matrix')
	return outfile

mds_r_program = """
d <- read.table("%s")
coords <- cmdscale(d, k=%d)
write.table(coords, file="%s", row.names=F, col.names=F)
"""

def mds(distmat, k, outfile=None):
	"""perform <k>-dimensional mds using distance matrix in <distmat>"""
	if not outfile: outfile = distmat + ".coord"
	if os.path.exists(outfile):
		debug('reusing MDS result: %s' % outfile)
		return outfile
	mds_p = mds_r_program % (distmat, k, outfile)
	n = NTF()
	n.write("%s" % mds_p)
	n.flush()
	cmd = "R CMD BATCH %s %s.Rout" % (n.name, n.name)
	os_run(cmd, msg='cannot run MDS')
	n.close()
	os.unlink(n.name + ".Rout")
	return outfile

def distances(db, ref_db, outfile=None):
	"""
	for each database compound, use search program to get the 
	distance from it to all reference compounds.
	"""
	if not outfile: outfile = ref_db + ".distances"
	if os.path.exists(outfile):
		debug('reusing big distance table file: %s' % outfile)
		return outfile
	cmd = "%s %s %s %s > %s" % (DB2DB_DISTANCE,  CDB, db, ref_db, outfile)
	os_run(cmd, msg='cannot run %s to get "puzzle body"' % DB2DB_DISTANCE)
	return outfile
	
def wait_qsub(job_id):
	"""wait for qsub to finish"""
	import time
	cmd = "qstat | grep %s > /dev/null" % job_id
	while os.system(cmd) == 0:
		sys.stderr.write("qsub jobs still running. --%s  \r" % time.ctime())
		time.sleep(60)
	sys.stderr.write("\n")

def line_select(input, output, base=1):
	"""select lines corresponding to sample queries"""
	queries = set([int(i) for i in file(TEST_QUERIES)])
	of = file(output, 'w')
	for i, line in enumerate(file(input)):
		if i + base in queries:
			of.write(line)
	return output

def merge_result(jobs, inputs, output):
	"""check and merge qsub job result"""
	total_t = 0
	failed = []
	try:
		os.unlink(output)
	except: pass
	for job in jobs:
		e = glob.glob(job + '.e*')
		if len(e) != 1:
			critical("number of STDERR file for %s is not 1" % job)
		try:
			(l, t, r) = file(e[0]).readlines()[-1].split()
			assert l == 'Time:'
		except:
			failed.append(job)
			info("output of %s is unexpected. Please retry this job" % job)
		else:
			total_t += float(t)

	if not failed:
		for input in inputs:
			os_run('cat %s.out >> %s' % (input, output))

		info("total time in solving puzzle: %s" % total_t)
		os_run("echo %s >> timing" % total_t)

	return failed

def binarize_coord(input, output, dim):
	"""take character-based coordinate file, and convert to binary format based
	on lshkit's matrix format"""
	cmd = "%s %s %s %d" % (COORD_TO_BINARY, input, output, dim)
	os_run(cmd)
	return output
	
def eucsearch(record1, record2, output):
	"""perform euclidean space search for evaluation purpose to see how
	embedding works"""
	cmd = "%s %s %s %d | gzip > %s" % (EUCSEARCHTOOL, record1, record2,
			MAX_EUCSEARCH_RESULTS, output)
	os_run(cmd)
	return output

def accuracy(n, k):
	"""perform evaluation to see how embedding works"""
	cmd = "%s %s eucsearch.%d-%d recall" % (EVALUATOR,
			CHEMICAL_SEARCH_RESULTS, n, k)
	os_run(cmd)

def indexed_search(record, output="indexed.gz",
		evaluation_out="indexed.performance"):
	"""perform indexed search"""
	info("running indexed search")
	cmd = INDEXED_SEARCH % record + " | gzip > " + output
	info("cmd: "+cmd)
	ret, stdout, stderr = os_run(cmd, stderr=STORED) 
	elapse = 0
	prompt = ">>Query time:"
	for i in stderr.splitlines():
		print i
		if i.startswith(prompt):
			elapse += float(i[len(prompt):].strip())
	info( "total time: %.1f" % elapse )
	cmd = "echo %s > index.search.timing" % elapse
	os_run(cmd)

	cmd = "%s %s %s > %s" % (INDEXED_SEARCH_EVALUATOR, 
			CHEMICAL_SEARCH_RESULTS, output, evaluation_out)
	os_run(cmd)
	return output

def main(n, k, per_file=20000, input=None, post_action=None, coord_ready=False):
	"""take an input reference sdf, and get the coordinates
	   n = sample size
	   k = dimensionality
	"""
	if not coord_ready:
		if not input:
			# generate the sub-database for reference compounds
			try: import hashlib as md5
			except: import md5
			import time
			prefix = md5.md5(str(time.time())).hexdigest()
			input = subset_db(n, outfile="%s.cdb" % prefix)

		# perform mds to get coordinate for reference compounds
		mds_file = mds(dist_mat(input), k)

		# second half of the puzzle
		puzzle_file = distances(IDDB, input)

		# build the puzzle file
		pz = file(puzzle_file)
		cntr = 0
		jobs = []
		inputs = []
		job_tmpl = 'job-%s' % input[:3]
		prev_f = None
		for i, line in enumerate(pz):
			if i % per_file == 0:
				cntr += 1
				of = "%s-%s-%s" % (n, k, cntr)
				if prev_f: prev_f.close()
				f = file(of, 'w')
				f.write("%d %d\n" % (k, n))
				f.close()
				cmd = "cat %s >> %s" % (mds_file, of)
				os_run(cmd)
				job = '%s-%s.sh' % (job_tmpl, cntr)
				cmd = "echo 'cd %s;%s %s' > %s" % (
					os.path.abspath(os.path.curdir), COORDTOOL, of, job)
				jobs.append(job)
				inputs.append(of)
				os_run(cmd)
				f = file(of, 'a')
				prev_f = f
			f.write(line)
		f.close()

		# use <post_action> to process jobs
		if post_action:
			jobs_todo = jobs
			for job in jobs_todo:
				if post_action == "dry":
					info("dry: job=" + job)
					continue
				elif post_action == "bash":
					os_run(post_action + " " + job+" 2> "+job+".e")
				else:
					os_run(post_action + " " + job)
			if post_action == 'qsub':
				wait_qsub(job_tmpl)
			# merge results
			coord_file = "coord.%s-%s" % (n, k)
			info("merging coordinate results")
			jobs_todo = merge_result(jobs, inputs, coord_file)
			if jobs_todo:
				print jobs_todo
				return

			coord_file2= "coord.query.%s-%s" % (n, k)
			line_select(coord_file, coord_file2)

	# making binary format of the coordinate files
	info("making binary version of coordinate files...")
	coord_file = "coord.%s-%s" % (n, k)
	coord_file2= "coord.query.%s-%s" % (n, k)
	r1 = binarize_coord(coord_file, "matrix.%s-%s" % (n,k), k)
	r2 = binarize_coord(coord_file2, "matrix.query.%s-%s" % (n,k), k)

	# perform euclidean search
	eucsearch(r1, r2, "eucsearch.%s-%s" % (n, k))
#	if input:
#		clean(input[:-3], n, k)

	# do accuracy test
	accuracy(n, k)

def query(r,d,query_sdf):
	t = time()
	assert os.path.isfile(query_sdf)
	work_dir = 'run-%s-%s' % (r, d)
	query_base = os.path.splitext(query_sdf)[0]
	query_cdb = os.path.join(work_dir,query_base+".cdb")

	info("creating query cdb")
	try:
		#produces dist.names?
		create_db(query_sdf,query_cdb,log_names=False,first=True)
		assert os.stat(query_cdb)[ST_SIZE] != 17
	except:
		print_exc()
		return False
	parsing_time = time() - t

	#create distance matrix
	info("creating distance matrix")
	query_dist = os.path.join(work_dir,query_base+".dist")
	try:
		unlink(query_dist)
	except:
		pass
	t=time()
	failed=False
	try:
		#uses ref.cdb
		comparer = DBComparer(CDB)
		comparer.compare(query_cdb,query_dist)
	except:
		failed=True
	if failed or os.stat(query_dist)[ST_SIZE] == 0:
		print_exc()
		error("Error in comparing to references")
		return False

	# solve the coordinate
	info("solving coordinate puzzle")
	query_coord = os.path.join(work_dir,query_base+".coord")
	try:
		#uses puzzle
		puzzle_file = os.path.join(work_dir,"puzzle")
		f = file(puzzle_file,'w')
		f.write('%s %s'%(d,r))
		f.close()
		os_run("cat %s >> %s" % (os.path.join(work_dir,'coord.%s-%s' % (r,d)),puzzle_file))
		solver = CoordinateSolver(puzzle_file)
		solverResult = solver.solve(query_dist).strip()
		assert solverResult
	except:
		print_exc()
		error("Error in embedding")
		return False
	f = file(query_coord,'w')
	f.write(solverResult)
	f.close()
	coord_time = time() - t

	# perform lsh search
	info("performing lsh search")
	t=time()
	try:
		#uses matrix
		lshsearcher = LSHSearcher(os.path.join(work_dir,"matrix.%s-%s" % (r,d)),lsh_param)
		lshResult = lshsearcher.search(solverResult).strip()
		assert lshResult
		print lshResult
	except:
		print_exc()
		error("Error in performing LSH query")
		info("solver result: \n"+solverResult)
		return False
	lsh_time = time() - t

	#refine
	info("refining")
	t = time()
	try:
		query_lsh = "candidates.data"
		f = file(query_lsh,'w')
		f.write(lshResult)
		f.close()

		#refiner = Refiner("db.fp_cdb",200)
		refiner = Refiner(CDB,d)
		os_run("cp "+query_cdb+" query.fp_cdb")	
		refineResult = refiner.refine(query_cdb+" "+query_lsh)
		assert refineResult
	except:
		print_exc()
		error("Error in performing refinement")
		return False
	refine_time = time() - t

	info("outputing results")
	sys.stderr.write('timing: parsing=%s embedding=%s lsh=%s refine=%s \n' %
		(parsing_time, coord_time, lsh_time, refine_time))

	names = [i.strip() for i in file(CDB+".names")]

	f = file(query_base+".out",'w')
	f.write('# %s %s %s %s\n' % 
		(parsing_time, coord_time, lsh_time, refine_time))
	for pair in refineResult.split():
		seq_id, dist = pair.split(':')
		cid = names[int(seq_id)-1]
		f.write('%s %s\n' %(cid,dist))
	f.close()

	info("done")
	return True



def time_and_result(line):
	if line.startswith('/t:'):
		time_block,line = line.split(None,1)
		time= float(time_block.split(':')[1])
	else:
		time = 0
	return (time,line)

def clean(input, n, k):
	info("cleaning up files")
	for i in glob.glob('%s.*' % input):
		info("removing %s" % i)
		os.unlink(i)
	if input:
		for i in glob.glob('job-%s-*' % input[:3]):
			info("removing %s" % i)
			os.unlink(i)
	for i in glob.glob('%s-%s-*' % (n, k)):
		info("removing %s" % i)
		os.unlink(i)
		

if __name__ == '__main__':
	per_file = 20000
	import optparse
	p = optparse.OptionParser(usage=
		"usage: %prog [option] [[clean] sdfname]\n"
		"       %prog [option] accuracy\n"
		"       %prog [option] accuracy++\n"
		" accuracy: perform accuracy analysis only using existing data\n"
		" accuracy++: pick up from the point coordiante calculation is done"
	)
	p.add_option("-r", help="number of references", dest="r")
	p.add_option("-d", help="number of dimensions", dest="d")
	p.add_option("-m", help="similarity measure to use", dest="m")
	p.add_option("-q", help="query", dest="q")
	p.add_option("--dry-run", help="dry run", dest="dry", action="store_true",
		default=False)
	p.add_option("-s", "--slice", help="number of puzzles per job", dest="s")
	opts, args = p.parse_args()

	if opts.m is not None and opts.m:
		DB2DB_DISTANCE = os.path.join(BINDIR, "%s.%s" % (DB2DB_DISTANCE,opts.m))


	if opts.r is None or opts.d is None:
		sys.stderr.write("must specify r and d. Use --help to see usage\n")
		sys.exit(1)
	try:
		r = int(opts.r)
		d = int(opts.d)
		if opts.s: per_file = int(opts.s) 
	except ValueError:
		sys.stderr.write("r, d and s must all be integer.\n")
		sys.exit(1)
	if opts.dry: post_action = "dry"
	else: post_action = processor


	if opts.q:
		query(r,d,opts.q)
	elif len(args) == 0:
		work_dir = 'run-%s-%s' % (r, d)
		os_run("mkdir -p %s" % work_dir)
		os.chdir(work_dir)
		main(r, d, per_file=per_file, post_action=post_action)
	elif len(args) == 1 and args[0] == "accuracy":
		# use accuracy analysis tool to analyze accuracy used existing data
		work_dir = 'run-%s-%s' % (r, d)
		if work_dir != os.path.basename(os.path.curdir):
			os.chdir(work_dir)
		accuracy(r, d)
	elif len(args) == 1 and args[0] == "accuracy++":
		# pick up from the point coordiate calculation is done
		work_dir = 'run-%s-%s' % (r, d)
		if work_dir != os.path.basename(os.path.curdir):
			os.chdir(work_dir)
		main(r, d, coord_ready=True)

	elif len(args) == 1 and args[0] == "indexed":
		# pick up from the point coordiate calculation is done
		work_dir = 'run-%s-%s' % (r, d)
		if work_dir != os.path.basename(os.path.curdir):
			os.chdir(work_dir)
		indexed_search('matrix.%d-%d' % (r, d))

	elif len(args) == 1:
		# use existing sub-database
		work_dir = os.path.dirname(args[0])
		input = os.path.basename(args[0])
		if work_dir: os.chdir(work_dir)
		main(r, d, input=input, per_file=per_file, post_action=post_action)
	else:
		# comamnd is "clean some.sdf"
		assert args[0] == 'clean'
		work_dir = os.path.dirname(args[1])
		input = os.path.basename(args[1])
		if work_dir: os.chdir(work_dir)
		clean(input, r, d)
