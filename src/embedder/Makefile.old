OBJS = solver.o routines.o ../atompair/profiling.o
EXEC = coord icoord coord_server
#CFLAGS = -g -D_TIMER_
#CXXFLAGS = -g -D_TIMER_
CFLAGS = -mtune=nocona  -O3  -fomit-frame-pointer -funroll-loops -mfpmath=sse  -D_TIMER_ -I ../atompair/
CXXFLAGS = -mtune=nocona -O3  -fomit-frame-pointer -funroll-loops -mfpmath=sse  -D_TIMER_ -I ../atompair/
LDFLAGS = -static

all : $(EXEC)
.c.o:
.cc.o:

coord : $(OBJS) coord.o
	$(CXX) -o $@ coord.o $(OBJS) -lm -lf2c -static
icoord : $(OBJS) icoord.o
	$(CXX) -o $@ icoord.o $(OBJS) -lm -lf2c -lreadline
coord_server : $(OBJS) coord_server.o
	$(CXX) -o $@ coord_server.o $(OBJS) -lm -lf2c 
	
clean :
	rm -f $(OBJS) $(EXEC)
