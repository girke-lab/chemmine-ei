import os
BASEDIR = os.path.dirname(os.path.abspath(__file__))

# path to the chemical database
chem_db = os.path.join(BASEDIR, 'data', 'chemdb.cdb')

# number of sample queries to create if not present
n_sample_queries = 1000		

# processor through what? use 'qsub' if you have access to cluster or use
# 'bash' for local processing
processor = 'qsub'

# parameter for LSH search. See 
# http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html
lsh_param = " -W 1.39564 -M 19 -L 30 -K 600 -S 30 -T 30 "
