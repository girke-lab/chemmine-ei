"""On-Demand Embedding. After the embedding of a whole database, you can embed
other compounds on a per-compound basis. This program is for this embedding.
You must specificy which reference compound database and reference coordinates
to use."""
import os
import sys
from ei import BINDIR, os_run, CDB
from tempfile import mkdtemp

def gen_subdb(ref_db_path, measure):
    """generate real database for iddb-style subdatabase"""
    ref_real_db = ref_db_path + '.db'
    from stat import ST_SIZE
    if os.path.isfile(ref_real_db) and os.stat(ref_real_db)[ST_SIZE]:
        sys.stderr.write("Reusing database " + ref_real_db)
    else:
        db_writer = os.path(BINDIR, 'db_subset')
        if measure: db_writer += ('.' + measure)
        os_run('%s %s %s' % (CDB, ref_db_path, ref_real_db), 
            msg="Cannot generate subdatabase")
    return ref_real_db

def embed(compound, r, d, ref_db_path, ref_coord, db_builder, db2db_distance):
    """embed <compound> in <d> dimensional space using <r> reference compounds
    specified in <ref_db_path>. db_builder is used to parse the input
    compound into a database file"""
    dir = mkdtemp()
    out_path = os.path.join(dir, os.path.basename(compound)) + '.cdb'
    os_run("%(cmd)s %(inp)s %(out)s" % dict(cmd=db_builder, inp=compound,
                                            out=out_path),
           msg="Cannot parse input file")

    # now build the puzzle file
    pz_path = os.path.join(dir, 'puzzle')
    f = file(pz_path, 'w')
    f.write("%d %d\n" % (d, r))
    f.flush()
    f.close()

    os_run("cat %s >> %s" % (ref_coord, pz_path))

    # compare the input to reference database
    os_run("%(cmd)s %(cmp)s %(ref_db)s >> %(out)s" % dict(cmd=db2db_distance,
                                                         cmp=out_path,
                                                         ref_db=ref_db_path,
                                                         out=pz_path),
           msg="cannot compare input to reference database")

    # solve the puzzle
    solver = os.path.join(BINDIR, 'coord')
    os_run("%(cmd)s %(inp)s" % dict(cmd=solver, inp=pz_path),
           msg="Cannot run embedder")
    f = file(pz_path + '.out')
    x = f.read()
    f.close()

    from shutil import rmtree
    rmtree(dir)
    return x

if __name__ == '__main__':
    # read command line arguments: <R> and <D>
    import optparse
    p = optparse.OptionParser(usage=
        "usage: %prog [option] compound.sdf")
    p.add_option('-m', help="similarity measure to use", dest="m")
    p.add_option('-x', help="reference compound database file", dest="ref")

    opts, args = p.parse_args()
    
    ref_db_path = opts.ref
    if not os.path.isfile(ref_db_path):
        sys.stderr.write("Invalid reference compound database: no such file.")
        sys.stderr.write("You specified\n")
        sys.stderr.write(ref_db_path)
        sys.exit(1)
    ref_coord = ref_db_path + '.distmat.coord'
    if not os.path.isfile(ref_coord):
        sys.stderr.write("Invalid reference coordinate file: no such file.")
        sys.stderr.write("Missing:\n")
        sys.stderr.write(ref_coord)
        sys.exit(1)
    
    f = file(ref_coord)
    d = len(f.readline().split())
    for i, _ in enumerate(f): pass
    r = 2 + i
    sys.stderr.write("r = %d d = %d\n" % (r, d))

    db_builder = os.path.join(BINDIR, 'db_builder')
    db2db_distance = os.path.join(BINDIR, 'db2db_distance')
    if opts.m:
        db_builder += ('.' + opts.m)
        db2db_distance += ('.' + opts.m)

    if not os.path.isfile(db_builder):
        sys.stderr.write("Cannot find database builder to parse your input")
        sys.stderr.write("\nI cannot find: ")
        sys.stderr.write(db_builder)
        sys.exit(1)

    ref_real_db = gen_subdb(ref_db_path, opts.m)
 
    print embed(args[0], r, d, ref_real_db, ref_coord, db_builder, db2db_distance)

# vim:tw=78:ts=4:sw=4:expandtab
