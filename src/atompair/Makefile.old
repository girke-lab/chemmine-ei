CXXFLAGS =      -g -Wall -fPIC -D_BUILD_DB_BUILD_ -DPRELOAD_DATABASE -O3 -mtune=nocona #-DDEBUG

OBJS =          molecule.o formats.o desc.o similarity.o simpledb.o search.o profiling.o 
FP_OBJS =          molecule.o formats.o desc.o simpledb.o search.o profiling.o
EXTRA_OBJS =    db_build.o compare.o db_search.o script.o db_search_range.o db_search_self.o db_isearch.o  db_compare_server.o db2db_distance.o

LIBS =    -lreadline -static

TARGETS =   libdescriptors.a  compare db_build db_search db_search_self db_search_range db2db_distance db_subset

all: $(TARGETS)

libdescriptors.a: $(OBJS)
	ar -r $@ $(OBJS)

compare: $(OBJS) compare.o 
	$(CXX) -o compare compare.o $(OBJS) $(LIBS)
db_build: $(OBJS) db_build.o
	$(CXX) -o db_build db_build.o $(OBJS) $(LIBS)

db_search: $(OBJS) db_search.o
	$(CXX) -o db_search db_search.o $(OBJS) $(LIBS)  
db_isearch: $(OBJS) db_isearch.o
	$(CXX) -o db_isearch db_isearch.o $(OBJS) $(LIBS)  

db_search_self: $(OBJS) db_search_self.o
	$(CXX) -o db_search_self db_search_self.o $(OBJS) $(LIBS)  

db_search_range: $(OBJS) db_search_range.o
	$(CXX) -o db_search_range db_search_range.o $(OBJS) $(LIBS)  
db2db_distance: $(OBJS) db2db_distance.o
	$(CXX) -o db2db_distance db2db_distance.o $(OBJS) $(LIBS)
db_subset: $(OBJS) db_subset.o
	$(CXX) -o $@ $@.o $(OBJS) $(LIBS)
        
clean:
	rm -f $(OBJS) $(EXTRA_OBJS) $(TARGETS)

dist : all
	rm -rf dist
	mkdir dist
	cp *.cc *.h Makefile Makefile.openbabel dist
	tar -czf descriptor.tar.gz dist
	rm -rf dist
