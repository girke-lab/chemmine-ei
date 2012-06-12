# Eddie's utility module

import sys
import os
from subprocess import Popen, PIPE, STDOUT
import logging
import logging.config
from tempfile import mkstemp, TemporaryFile, NamedTemporaryFile
from stat import ST_SIZE
from copy import copy

JOB_DIR_PREFIX = '~www/dockgui' # Where Job folders are located
JOB_DIR_PREFIX = os.path.expanduser(JOB_DIR_PREFIX)

p_dir = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(p_dir, "elog.config"))
elog = logging.getLogger('elog')
configFile = os.path.join(os.path.expanduser('~'),".eirc")

not_empty = lambda x: os.path.exists(x) and os.stat(x)[ST_SIZE]

class GzipFile(object):
    def __init__(self, path):
        self.p = Popen('zcat %s' % path, bufsize=1024, close_fds=True,
                shell=True, stdout=PIPE, stderr=None)
        self.status = 'open'
    
    def __iter__(self):
        assert self.status == 'open'
        for i in self.p.stdout:
            yield i

    def next(self):
        assert self.status == 'open'
        return self.p.stdout.next()

    def close(self):
        try:
            self.p.terminate()
        except:
            from signal import SIGTERM
            try:
                os.kill(self.p.pid, SIGTERM)
            except: pass
        self.p.stdout.close()
        self.status = 'closed'

def getConfig():
	"write out a default config file in ~/.eirc if not found, else, return the existing file"
	if(not os.path.exists(configFile)):
		out = open(configFile,"w")
		out.write("""
import os
BASEDIR =  os.path.abspath(".")

# path to the chemical database                    
chem_db = os.path.join(BASEDIR, 'data', 'chemdb.cdb')

# number of sample queries to create if not present
n_sample_queries = 1000

# processor through what? use 'qsub' if you have access to cluster or use
# 'bash' for local processing
processor = 'bash'

K=600 # number of nearest neighbors in embedded space

# parameter for LSH search. See
# http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html
lsh_param = " -W 1.39564 -M 19 -L 30 -K %s -T 30 " % K
		""")
		out.close()
	return configFile



def symlink_all(dir):
    """
    symlink all files and dirs(not recursive) under dir to the current dir
    """
    import glob
    for i in glob.glob(os.path.join(dir, '*')):
        os.symlink(i, os.path.basename(i))
    
def is_float(s):
    if not s: return False
    try:
        float(s)
    except ValueError:
        return False
    return True
def pdb_residue_iter(pdb):
    """
    interate residues in a pdb file. pdb is a file-like object
    """
    last_rid = None
    buf = ''
    for line in pdb:
        section = line[:6].strip()
        if section != 'ATOM' and section != 'HETATM':
            continue 
        if not is_float(line[30:38]) or not is_float(line[38:46]) \
                                                or not is_float(line[46:54]):
            continue

        rid = int(line[22:26])
        if rid != last_rid:
            if buf:
                yield {'':buf, 'chain':chain, 'resname':resname,
                                                                'rid':last_rid}
            buf = ''
        buf += line
        chain, resname = line[21], line[17:20] 
        last_rid = rid

    if buf:
        yield {'':buf, 'chain':chain, 'resname':resname, 'rid':last_rid}
 
def mol2_iter(input):
    """
    iterate a mol2 file. input is a file-like object
    """
    sep = r'@<TRIPOS>MOLECULE'
    assert input.next().strip() == sep
    buf = sep + '\n'

    for i in input:
        if i.strip() == sep:
            yield buf
            buf = sep + '\n'
        else:
            buf += i
    yield buf

def tmpfile(*args, **kargs):
    h, path = mkstemp(*args, **kargs)
    os.close(h)
    return path

total_time = dict()
def timed(func):
    def f(*args, **kargs):
        if '_timer_name_' not in kargs: _timer_name_ = func.__name__
        else: _timer_name_ = kargs.pop('_timer_name_')
        from time import time
        if _timer_name_ not in total_time: total_time[_timer_name_] = 0

        __start_time__ = time()
        ret = func(*args, **kargs)
        elapse = time() - __start_time__

        elog.debug("time on %s is %f" % (_timer_name_, elapse))
        total_time[_timer_name_] += elapse

        return ret
    return f

def open_argument_file(ext=''):
    """
    a decorator to allow a function to transparently take an argument of
    either a path or a file-like object. This argument must be the first
    argument. A __with_path argument will be added if accepted to point to
    the path of the file if available.
    """
    def deco(func):
        def f(f, *args, **kargs):
            if isinstance(f, str):
                if not f.endswith(ext):
                    raise Exception('Expect input to be a %s file' % ext)
                fo = file(f)
                kargs['__with_path'] = f
            else:
                fo = f
            try:
                return func(fo, *args, **kargs)
            except TypeError, e:
                if "got an unexpected keyword argument '__with_path'" in str(e):
                    del kargs['__with_path']
                    return func(fo, *args, **kargs)
                else:
                    raise
        return f
    return deco
        
def get_job_dir(job_id):
    """
    given the job_id, return the job dir
    """
    job_id_name = "%d" % int(job_id)
    return os.path.join(JOB_DIR_PREFIX, '0' + job_id_name[-1:], job_id_name)

DISCARDED = -1
DEFAULT = 0
STORED = 1
JOINED = 2

class ExternalProgramError(Exception):
    pass

class OS_Runner(object):
    def __init__(self, **kargs):
        self.defaults = kargs

    def update(self, **kargs):
        self.defaults.update(kargs)

    def __call__(self, cmd, **kargs):
        params = copy(self.defaults)
        params.update(kargs)
        return os_run2(cmd, **params)

def os_run2(cmd, msg=None, exception=True, prefix=None, stdout=DEFAULT,
            stderr=DEFAULT, output=sys.stdout):
    """
    new generation of os_run, that invokes a UNIX command, waits it to finish
    and returns the return code. The return is always a three-element tuple:
        (return_code, stdout_str, stderr_str)

    - msg is the message used in Exceptions and logs.
    - when exception is set, an Exception will be raised when return code is 
    non-zero.
    - prefix is used when stdout and stderr is printed.
    - stdout and stderr controls the handling of stdout and stderr. Possible 
    values are:
        - DEFAULT: no handling. print to terminal, subject to prefix
        - STROED: stored in a string and returned in the tuple
        - JOINED: for stderr only, sent to stdout
        - DISCARDED: discarded.
    """
    if prefix is None: prefix = ''
    elog.debug("invoking: %s" % cmd)
    null = file('/dev/null', 'w')
    if stdout == DEFAULT:
        _stdout = PIPE
    elif stdout == STORED:
        _stdout = TemporaryFile()
    elif stdout == DISCARDED:
        _stdout = null
    else:
        raise ExternalProgramError("Unknown mode '%s' passed to os_run" % stdout)

    if stderr == DEFAULT:
        _stderr = None
    elif stderr == JOINED:
        _stderr = STDOUT
    elif stderr == STORED:
        _stderr = TemporaryFile()
    elif stderr == DISCARDED:
        _stderr = null
    else:
        raise ExternalProgramError("Unknown mode '%s' passed to os_run" % stderr)

    try:
        subp = Popen(cmd, bufsize=1024, close_fds=True, shell=True,
                    stdout=_stdout, stderr=_stderr)

        if _stdout == PIPE:
            for i in subp.stdout:
                output.write(prefix + i)
        subp.wait()

        ret_out, ret_err = None, None
        if stdout == STORED:
            _stdout.seek(0)
            ret_out = _stdout.read()
            _stdout.close()
        if stderr == STORED:
            _stderr.seek(0)
            ret_err = _stderr.read()
            _stderr.close()

        if exception:
            if subp.returncode != 0:
               raise ExternalProgramError('Error in running %s, return code: %d' % (cmd, subp.returncode))

        return (subp.returncode, ret_out, ret_err)
    except OSError:
        if not msg: msg = 'ERROR: cannot run %s' % cmd
        critical(msg)
        raise ExternalProgramError

def os_run(cmd, msg=None, exit_on_fail=True, pipe_out=False,
        return_out=False, stdout_prefix="[default]",
        output=sys.stdout):
    """
    invoke a UNIX command. print `msg' on error.

    if return_out is set, output will be fully read and returned as a tuple
    (returncode, stdout)

    pipe_out is deprecated
        if pipe_out is set, a file-like object containing the stdout will be
        returned.  Only if pipe_out is False will this function return the
        returncode. As a hack, you can pass in a (non-empty) list, and
        pipe_out will put the Popen object in the first slot of the list.
        Dirty hack!

    Otherwise, the function returns returncode
    """
    if msg: elog.debug(msg)
    elog.debug("invoking: %s" % cmd)

    if stdout_prefix == '[default]':
        stdout_prefix = cmd.split(' ', 1)[0] + ":\t"

    try:
        subp = Popen(cmd, bufsize=1024, close_fds=True, shell=True,
                    stdout=PIPE)
        if return_out:
            stdout, _ = subp.communicate()
            return (subp.returncode, stdout)
        elif pipe_out:
            if isinstance(pipe_out, list): pipe_out[0] = subp
            return subp.stdout
        else:
            for i in subp.stdout:
                output.write("%s%s" % (stdout_prefix, i))
            subp.wait()
            return subp.returncode
    except OSError:
        if not msg: msg = 'ERROR: cannot run %s' % cmd
        critical(msg)
        if exit_on_fail: sys.exit(1)
        else: raise Exception(msg)

def run_or_raise(cmd, exception, *args, **kargs):
    if os_run(cmd, *args, **kargs) != 0:
        raise exception

def walk_dir(path_segs, dir=".", level=0):
    """
    build all possible paths based on the path segments in path_segs.
    This function returns an iterator that can be used in a for loop.
    
    example path_segs would be:
    path_segs = [
                    [
                        ('dir', 'run.ligand'),
                        ('file', 'rec.pdb'),
                        ('file', 'poc.txt'),
                    ],
                    [
                        ('dir', '1.A'),
                        ('dir', '1.B'),
                        ('dir', '2.A'),
                        ('dir', '2.B'),
                    ],
                    [
                        ('file', 'a.energy', 'hint1'),
                        ('file', 'a.new.mol2'),
                    ],
            ]
    """
    hint = None
    for path_info in path_segs[level]:
        _type, path_seg = path_info[:2]
        if len(path_info) == 3: hint = path_info[2]
        else: hint = None
        if _type == 'dir':
            yield (_type, os.path.join(dir, path_seg), level, hint)
            if level < len(path_segs) - 1:
                for sub in walk_dir(path_segs, os.path.join(dir, path_seg),
                                    level + 1):
                    yield sub
        else:
            yield (_type, os.path.join(dir, path_seg), level, hint)

## RMSD COMMANDS ##
convert_cmd = """bash -c "convert.py --no_hyd --i=%s --o=%s" """
convert_cmd2 = """ file2file.py -c %s %s """
#rmsd_cmd = """ fconv.linux %s -rmsd --s=%s """
rmsd_cmd = """ rmsd %s %s """

def canonical(structure):
    """
    canonically convert structure to a pdb file
    """
    g = tmpfile(suffix='.mol2')
    os_run2(convert_cmd % (structure, g), stdout=DISCARDED, stderr=DISCARDED)
    g2 = tmpfile(suffix='.pdb')
    os_run2(convert_cmd2 % (g, g2), stdout=DISCARDED, stderr=DISCARDED)
    os.unlink(g)
    return g2

def rmsd(s1, s2, skip_canon=False):
    if not skip_canon:
        s1 = canonical(s1)
        s2 = canonical(s2)
    if not s1.endswith('pdb') or not s2.endswith('pdb'):
        raise ExternalProgramError("RMSD calculation only works on PDB")
    ret, out, _ = os_run2(rmsd_cmd % (s1, s2), stdout=STORED)
    if not skip_canon:
        os.unlink(s1)
        os.unlink(s2)
    if 'Error!' in out:
        raise ExternalProgramError(out.splitlines()[-1])
    return float(out.splitlines()[-1].split()[-1])
#    rmsd = None
#    for line in out.splitlines():
#        _rmsd = float(line.split()[-1])
#        if rmsd is None or _rmsd < rmsd:
#            rmsd = _rmsd
#    return rmsd
## ENDS RMSD ##

@open_argument_file(ext='.mol2')
def fix_mol2_names(f, names, saveas=None, __with_path=None, inplace=False):
    """
    fix names of molecules in file mol2 using a list <names>
    It will do it in place if saveas is not set

    if <names> is not a list, it will apply it to all molecules
    """
    t = NamedTemporaryFile()
    cntr = 0
    for i, mol2 in enumerate(mol2_iter(f)):
        lines = mol2.splitlines()
        if isinstance(names, list):
            try:
                lines = [lines[0]] + [names[i]] + lines[2:]
            except IndexError:
                elog.error("Running out of names. Will only keep the finished %s molecules" % i)
                break
        else:
            lines = [lines[0]] + [names] + lines[2:]
        t.write('\n'.join(lines))
        t.write('\n')
        cntr = i + 1
    t.flush()

    if inplace:
        saveas = __with_path
    elif not saveas and __with_path:
        saveas = __with_path
        elog.warning("Warning: Overwriting input file")

    elog.info("%d molecules fixed" % cntr)
    if isinstance(names, list):
        if len(names) > cntr:
            elog.warning("Warning: There are unused names in the list")

    from shutil import copy
    copy(t.name, saveas)
    t.close()

import cStringIO,operator

def indent(rows, hasHeader=False, headerChar='-', delim=' | ', justify='left',
           separateRows=False, prefix='', postfix='', wrapfunc=lambda x:x):
    """Indents a table by column.
       - rows: A sequence of sequences of items, one sequence per row.
       - hasHeader: True if the first row consists of the columns' names.
       - headerChar: Character to be used for the row separator line
         (if hasHeader==True or separateRows==True).
       - delim: The column delimiter.
       - justify: Determines how are data justified in their column. 
         Valid values are 'left','right' and 'center'.
       - separateRows: True if rows are to be separated by a line
         of 'headerChar's.
       - prefix: A string prepended to each printed row.
       - postfix: A string appended to each printed row.
       - wrapfunc: A function f(text) for wrapping text; each element in
         the table is first wrapped by this function."""
    # closure for breaking logical rows to physical, using wrapfunc
    def rowWrapper(row):
        newRows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in map(None,*newRows)]
    # break each logical row into one or more physical ones
    logicalRows = [rowWrapper(row) for row in rows]
    # columns of physical rows
    columns = map(None,*reduce(operator.add,logicalRows))
    # get the maximum of each column by the string length of its items
    maxWidths = [max([len(str(item)) for item in column]) for column in columns]
    rowSeparator = headerChar * (len(prefix) + len(postfix) + sum(maxWidths) + \
                                 len(delim)*(len(maxWidths)-1))
    # select the appropriate justify method
    justify = {'center':str.center, 'right':str.rjust, 'left':str.ljust}[justify.lower()]
    output=cStringIO.StringIO()
    if separateRows: print >> output, rowSeparator
    for physicalRows in logicalRows:
        for row in physicalRows:
            print >> output, \
                prefix \
                + delim.join([justify(str(item),width) for (item,width) in zip(row,maxWidths)]) \
                + postfix
        if separateRows or hasHeader: print >> output, rowSeparator; hasHeader=False
    return output.getvalue()

# written by Mike Brown
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/148061
def wrap_onspace(text, width):
    """
    A word-wrap function that preserves existing line breaks
    and most spaces in the text. Expects that existing line
    breaks are posix newlines (\n).
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line[line.rfind('\n')+1:])
                         + len(word.split('\n',1)[0]
                              ) >= width)],
                   word),
                  text.split(' ')
                 )

def gen_subdb(ref_db_path, measure,db_writer,main_cdb):
    """generate real database for iddb-style subdatabase"""
    ref_real_db = ref_db_path + '.db'
    from stat import ST_SIZE
    if os.path.isfile(ref_real_db) and os.stat(ref_real_db)[ST_SIZE]:
        sys.stderr.write("Reusing database " + ref_real_db)
    else:
        if measure: db_writer += ('.' + measure)
        os_run('%s %s %s %s' % (db_writer,main_cdb, ref_db_path, ref_real_db), 
            msg="Cannot generate subdatabase")
    return ref_real_db


import re
def wrap_onspace_strict(text, width):
    """Similar to wrap_onspace, but enforces the width constraint:
       words longer than width are split."""
    wordRegex = re.compile(r'\S{'+str(width)+r',}')
    return wrap_onspace(wordRegex.sub(lambda m: wrap_always(m.group(),width),text),width)

import math
def wrap_always(text, width):
    """A simple word-wrap function that wraps text on exactly width characters.
       It doesn't split the text in words."""
    return '\n'.join([ text[width*i:width*(i+1)] \
                       for i in xrange(int(math.ceil(1.*len(text)/width))) ])
    

# vim:expandtab
