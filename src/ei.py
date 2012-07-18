#!/usr/bin/env python
from logging import info, warning, error, debug, critical, root, NOTSET
import sys
import os
import random
from tempfile import NamedTemporaryFile as NTF,mkdtemp
import glob
from io import StringIO
from subprocess import Popen, PIPE, STDOUT,check_call
import re
from time import time
from stat import ST_SIZE
from traceback import print_exc
from eutils import OS_Runner, STORED, DISCARDED, getConfig,gen_subdb, time_function
from eutils.sdfiterator import sdf_iter
from eutils.coord import CoordinateSolver
from eutils.fpdbcompare import DBComparer
from eutils.lshsearch import LSHSearcher
from eutils.refineserver import Refiner


root.setLevel(NOTSET)
os_run = OS_Runner()

BINDIR = "" 
BASEDIR = "."
DATADIR = "data"
DB2DB_DISTANCE = os.path.join(BINDIR, "ei-db2db_distance")
DB_SUBSET = os.path.join(BINDIR,"ei-db_subset")
DB_BUILDER = os.path.join(BINDIR,"ei-db_builder")
EUCSEARCHTOOL = os.path.join(BINDIR, "ei-euclid_search")
EVALUATOR = os.path.join(BINDIR, "ei-evaluator")
COORD_TO_BINARY = os.path.join(BINDIR, "ei-bin_formatter")
INDEXED_SEARCH_EVALUATOR = os.path.join(BINDIR, "ei-comparesearch")
COORDTOOL = os.path.join(BINDIR, "ei-coord")
INDEXED_SEARCH = os.path.join(BINDIR, "ei-isearch")
SINGLE_SEARCH = os.path.join(BINDIR,"ei-single_search")
DB_ISEARCH = os.path.join(BINDIR,"ei-db_isearch")
K = 600

execfile(getConfig())

BASEDIR =  os.path.abspath(BASEDIR)
DATADIR = os.path.join(BASEDIR, DATADIR)

localConfig = os.path.join(DATADIR,"eirc")
if os.path.exists(localConfig):
	execfile(localConfig)

CDB = os.path.join(DATADIR, 'chem.db')
IDDB = os.path.join(DATADIR, 'main.iddb')
TEST_QUERIES = os.path.join(DATADIR, 'test_query.iddb')
CHEMICAL_SEARCH_RESULTS = os.path.join(DATADIR, 'chemical-search.results.gz')
MAX_EUCSEARCH_RESULTS = 50000

lsh_param = "%s -K %d" % (lsh_param,K)
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
#		error("!!!!! PLEASE GENERATE THE CHEMICAL SEARCH RESULT !!!!!")
#		error("!!!!! EXAMPLE:                                   !!!!!")
#		error("      ei-db_search -id chem.db main.iddb test_query.iddb 50000 | gzip > chemical-search.results.gz")
#		error(" ")
#		error("Please press any key to continue or CTRL-C to exit")
#		sys.stdin.read()
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
	print("dist cmd: "+cmd)
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

def gen_chemical_search_results():
	
	from gzip import GzipFile as zfile
	if not os.path.isfile(CHEMICAL_SEARCH_RESULTS):
		f = zfile(CHEMICAL_SEARCH_RESULTS,"w")
		distances = runDb2Db(
				  os.path.join(DATADIR,"chem.db"),
				  os.path.join(DATADIR,"test_query.iddb"),
				  os.path.join(DATADIR,"main.iddb"))
		f.write("\n".join(
			[" ".join(
				[ "%d:%f" % pair for pair in
					bestCandidates(queryDist,None,50000) ]) 
			for queryDist in distances]))
		f.close()
				  
		#check_call("ei-db_search -id %s %s %s  50000 | gzip > %s" % 
				#( os.path.join(DATADIR,"chem.db"),
				  #os.path.join(DATADIR,"main.iddb"),
				  #os.path.join(DATADIR,"test_query.iddb"),
				  #CHEMICAL_SEARCH_RESULTS),
				#shell=True)

def accuracy(n, k):
	"""perform evaluation to see how embedding works"""
	gen_chemical_search_results()
	eucsearch("matrix.%s-%s"%(n,k),"matrix.query.%s-%s"%(n,k), "eucsearch.%s-%s" % (n, k))
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

	gen_chemical_search_results()
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

# don't test by default
#	# perform euclidean search
#	eucsearch(r1, r2, "eucsearch.%s-%s" % (n, k))

	if input:
		clean(input[:-3], n, k)

# don't test by default
#	# do accuracy test
#	accuracy(n, k)

def createQueryCdb(query_sdf,query_cdb):
	t = time()

	subp = Popen([DB_BUILDER,query_sdf,query_cdb], stderr=PIPE)
	num_compounds=-1
	for line in subp.stderr:
		m = re.match(r"(\d+) compounds read.",line)
		if m:
			num_compounds=int(m.group(1))
	if num_compounds == -1:
		warning("could not read how many compounds we have")

	return (time() - t,num_compounds)

def createDistanceMatrixBatch(ref_db,query_cdb,query_dist):
	t = time()

	comparer = DBComparer(ref_db)
	comparer.compare(query_cdb,query_dist)
	assert os.stat(query_dist)[ST_SIZE] != 0
	return time() - t

def solvePuzzle(r,d,ref_db,query_cdb,coord_file,puzzle_file):

	createPuzzleFile(r,d,coord_file,puzzle_file)

	os_run("%s %s %s >> %s" % (DB2DB_DISTANCE, query_cdb, ref_db, puzzle_file),
           msg="cannot compare input to reference database")

	os_run("%(cmd)s %(inp)s" % dict(cmd=COORDTOOL, inp=puzzle_file), msg="Cannot run embedder")
	if os.path.getsize(puzzle_file+".out") == 0:
		warning("empty "+puzzle_file+".out")
	f = file(puzzle_file + '.out')
	x = f.read()
	f.close()
	info("puzzle: "+x)
	return x

def createPuzzleFile(r,d,coord_file,puzzle_file):
	f = file(puzzle_file,'w')
	f.write('%s %s\n'%(d,r))
	f.close()

	os_run("cat %s >> %s" % (coord_file,puzzle_file))

def solvePuzzleBatch(r,d,query_dist,coord_file,puzzle_file):
	t = time()
	createPuzzleFile(r,d,coord_file,puzzle_file)
	solver = CoordinateSolver(puzzle_file)
	solverResult = solver.solve(query_dist).strip()
	assert solverResult

	return (time() - t,solverResult)

def lshSearch(matrix_file,solverResult):
	if not solverResult:
		warning("no solver result")
	f = file("coords.in","w")
	f.write(solverResult)
	f.close()

	subp = Popen("%s %s -D %s -C coords.in " % (SINGLE_SEARCH, lsh_param,matrix_file), shell=True, stdout=PIPE)

	return okResult(subp.stdout.read())
def batchQuery(outf,r,d,ref_db,queries,coord_file,matrix_file,names ):

	query_dist="query.dist"
	query_cdb="query.cdb"
	query_sdf="query.sdf"
	query_name_file="query.cdb.names"
	puzzle_file="puzzle"
	candidate_file="candidates.data"

	createPuzzleFile(r,d,coord_file,puzzle_file)

	parsing_time = 0
	dist_time,comparer = time_function(DBComparer,ref_db)
	coord_time,solver = time_function(CoordinateSolver,puzzle_file)
	lsh_time,lshsearcher = time_function(LSHSearcher,matrix_file,lsh_param)
	refine_time,refiner = time_function(Refiner,CDB,K)

	debug("names: "+str(names))
	for i,current_query in enumerate(sdf_iter(queries)):
		#write current query to file and convert to a cdb. will create a names file as well
		f = file(query_sdf,'w')
		f.write(current_query)
		f.close()
		parsing_time += time_function(createQueryCdb,query_sdf,query_cdb)[0]

		#read the name of this query
		f=file(query_name_file,'r')
		name=f.readline().rstrip()
		f.close()
		info("=============== "+name+" =====================")

		#perform search
		dist_time += time_function(comparer.compare,query_cdb,query_dist)[0]
		t,solverResult = time_function(solver.solve,query_dist)
		coord_time += t
		solverResult = solverResult.strip()

		t,lshResult = time_function(lshsearcher.search,solverResult)
		lsh_time += t
		lshResult = lshResult.strip()

		info("lshResult %d: %s" % (i,lshResult))

		distances = refineLocal(query_cdb,lshResult)

		for candidate_index,dist in distances:
			outf.write("%s\t%s\t%s\n" %(name,names[candidate_index],dist))
	
	sys.stderr.write('timing: parsing=%s embedding=%s lsh=%s refine=%s \n' %
				(parsing_time, dist_time + coord_time, lsh_time, refine_time))

def lshSearchBatch(matrix_file,solverResult):
	t = time()

	lshsearcher = LSHSearcher(matrix_file,lsh_param)
	lshResult = lshsearcher.search(solverResult).strip()
	assert lshResult

	f = file("candidates.data",'w')
	f.write(lshResult)
	f.close()
	return time() - t

def distToCandidates(candidate_indcies,query_db):
	"""return a matrix of distances, queries are rows, candidates are columns"""
	if candidate_indcies == None:
		debug( "writing candidates: "+str(candidate_indcies))
		target_file="candidates.db"
		f = file("candidates.iddb","w")
	
		for index in sorted(candidate_indcies):
			f.write("%d\n" % (index))
		f.close()

		info("call db_subset")
		check_call([DB_SUBSET,CDB,"candidates.iddb",target_db])
	else:
		target_db = CDB
	return runDb2Db(query_db,target_db)

def runDb2Db(*args):
	info("call db2db_distance")
	subp = Popen( [DB2DB_DISTANCE]+list(args),stdout=PIPE)
	t=[ [float(value) for value in line.strip().split()] for line in subp.stdout.read().splitlines()]
	#debug("done reading: "+str(t))
	return t

def bestCandidates(distances,candidate_indcies,num_neighbors=K):
	"""takes a one dimentional distance array and candidates indcies. 
		returns (candidate_index, distance) tuples"""
	#debug( "distances: "+str(distances))
	if candidate_indcies == None:
		candidate_indcies = range(1,cdbsize()+1)
		#debug("candidates: "+str(candidate_indcies))
		
	distances = [ (candidate_indcies[di[0]],di[1]) for di in enumerate(distances)]
	distances.sort(key=lambda x:x[1] )
	return distances[:num_neighbors]

def refineLocal(query_cdb,candidates):
	#debug("orig candidates: "+candidates)
	candidate_indcies = [int(s.split(":")[0]) for s in candidates.split() ]
	if not candidate_indcies:
		return []
	return bestCandidates( distToCandidates(candidate_indcies,query_cdb)[0], candidate_indcies)

def refine(query_cdb,candidates):
	f = file("candidates.data",'w')
	f.write(candidates)
	f.close()
	info("candidates: "+candidates)
	subp = Popen( [DB_ISEARCH,CDB,str(K)],stdin=PIPE,stdout=PIPE)
	subp.stdin.write(query_cdb+" candidates.data\n")
	subp.stdin.close()
	return okResult(subp.stdout.read())

def okResult(output):
	#info("output: "+output)
	result=None
	m = re.search(r"OK:(.*)\n",output)
	if m:
		return m.group(1)
	else:
		raise StandardError("no results found. output: "+output)

def refineBatch(query_cdb):
	t = time()

	refiner = Refiner(CDB,K)
	#os_run("cp "+query_cdb+" query.fp_cdb")	
	refineResult = refiner.refine("%s %s" %(query_cdb,"candidates.data"));
	assert refineResult

	return (time() - t,refineResult)

def query(r,d,query_sdf,ref_iddb):

	current_dir = os.path.abspath(".")
	query_sdf = os.path.join(current_dir,query_sdf)
	ref_iddb = os.path.join(current_dir,ref_iddb)
	work_dir = os.path.join(current_dir,'run-%s-%s' % (r, d))

	if not os.path.isfile(query_sdf):
		raise StandardError("query file "+query_sdf+" not found")
	if not os.path.isfile(ref_iddb):
		raise StandardError("reference file "+ref_iddb+" not found")
	if not os.path.isdir(work_dir):
		raise StandardError("working directory "+work_dir+" not found")

	temp_dir=mkdtemp()
	query_base = os.path.splitext(os.path.basename(query_sdf))[0]

	query_cdb = os.path.join(temp_dir,query_base+".cdb")

	matrix_file = os.path.join(work_dir,"matrix.%s-%s" % (r,d))
	coord_file = ref_iddb+".distmat.coord"

	try:
		parsing_time,num_compounds = createQueryCdb(query_sdf,query_cdb)
		info("found %s compounds" % num_compounds)
		ref_db = gen_subdb(ref_iddb,None,DB_SUBSET,CDB)
		names = [i.strip() for i in file(CDB+".names")]

		f = file(query_base+".out",'w')

		os.chdir(temp_dir)
		if num_compounds > 1:
			batchQuery(f,r,d,ref_db,query_sdf,coord_file,matrix_file,names)
			

			
			#dist_time = createDistanceMatrixBatch(ref_db,query_cdb,query_dist)
			#coord_time,solverResult = solvePuzzleBatch(r,d,query_dist,coord_file,puzzle_file)
			#lsh_time = lshSearchBatch(matrix_file,solverResult)
			#refine_time,refineResult = refineBatch(query_cdb)
			#sys.stderr.write('timing: parsing=%s embedding=%s lsh=%s refine=%s \n' %
				#(parsing_time, dist_time + coord_time, lsh_time, refine_time))
			#f.write('# %s %s %s %s\n' % 
				#(parsing_time, coord_time, lsh_time, refine_time))
		else:
			query_dist = os.path.join(temp_dir,query_base+".dist")
			puzzle_file = os.path.join(temp_dir,"puzzle")

			#refineResult = refine(query_cdb,lshSearch(matrix_file,solvePuzzle(r,d,ref_db,query_cdb,coord_file,puzzle_file)))
			#for pair in refineResult.split():
				#seq_id, dist = pair.split(':')
				#cid = names[int(seq_id)-1]
				#f.write('%s %s\n' %(cid,dist))

			refineResult = refineLocal(query_cdb,lshSearch(matrix_file,solvePuzzle(r,d,ref_db,query_cdb,coord_file,puzzle_file)))
			for seq_id,dist in refineResult:
				f.write('%s %s\n' %(names[int(seq_id)-1],dist))



			info("refine result: "+str(refineResult))

		f.close()
	
		os.chdir(current_dir)

		from shutil import rmtree
		#warning("NOT CLEANING UP")
		rmtree(temp_dir)
		
		
	except:
		print_exc()

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
		
def init(input_db,m):
	if not os.path.isdir("data"):
		os.mkdir("data")
	if m is not None and m:
		f=file(localConfig,"w")
		f.write("DB2DB_DISTANCE = '"+DB2DB_DISTANCE+"'\n")
		f.write("DB_SUBSET = '"+DB_SUBSET+"'\n")
		f.write("DB_BUILDER = '"+DB_BUILDER+"'\n")
		f.close()

	chemdb=os.path.join("data","chem.db")	
	if not os.path.isfile(chemdb):
		t,numCompounds = createQueryCdb(input_db,chemdb)
		f=file(os.path.join("data","main.iddb"),"w")
		for i in range(1,numCompounds+1):
			f.write(str(i)+"\n")
		f.close()

def check_data_exists():
	if not os.path.isfile(CDB):
		sys.stderr.write("No chem.db file found. Please create a database with the --init option\n")
		sys.exit(1)

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
	p.add_option("-x", help="reference cdb file", dest="x")
	p.add_option("--dry-run", help="dry run", dest="dry", action="store_true",
		default=False)
	p.add_option("-s", "--slice", help="number of puzzles per job", dest="s")
	p.add_option("--init", help="create initial database", dest="initdb")
	opts, args = p.parse_args()

	if opts.m is not None and opts.m:
		DB2DB_DISTANCE = os.path.join(BINDIR, "%s.%s" % (DB2DB_DISTANCE,opts.m))
		DB_SUBSET = os.path.join(BINDIR,"%s.%s" % (DB_SUBSET,opts.m))
		DB_BUILDER = os.path.join(BINDIR,"%s.%s" % (DB_BUILDER,opts.m))


	if opts.initdb is not None:
		init(opts.initdb,opts.m)
		sys.exit(0)

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


	check_data_exists()

	if opts.q and opts.x:
		query(r,d,opts.q,opts.x)
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
