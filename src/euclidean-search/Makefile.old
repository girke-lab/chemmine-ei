OBJS = recordfile.o search.o ../atompair/desc.o ../atompair/formats.o ../atompair/molecule.o ../atompair/similarity.o ../atompair/simpledb.o  ../atompair/profiling.o
CXXFLAGS = -O3 -mtune=nocona -D_REFINE_ -I ../atompair #-g
LDFLAGS = -static

EXEC = rnn knn scan self knn-self split
all: $(EXEC)

rnn : $(OBJS) rnn.o
	$(CXX) -o $@ $(OBJS) $@.o $(LDFLAGS)

knn : $(OBJS) knn.o
	$(CXX) -o $@ $(OBJS) $@.o $(LDFLAGS) 

knn-self : $(OBJS) knn-self.o
	$(CXX) -o $@ $(OBJS) $@.o $(LDFLAGS)

scan : $(OBJS) scan.o
	$(CXX) -o $@ $(OBJS) $@.o $(LDFLAGS)

self : $(OBJS) self.o
	$(CXX) -o $@ $(OBJS) $@.o $(LDFLAGS)
split : split.o
	$(CXX) -o $@ $@.o $(LDFLAGS)

clean :
	rm -f $(OBJS) $(EXEC) rnn.o knn.o knn-self.o scan.o self.o
