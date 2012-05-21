include ../Makefile.inc

CC=$(CXX)

# change this to help me locate header files for boost and boost-threadpool
#CXXFLAGS += -I../include -I ../../atompair/ -I/home/ycao/opt/src/boost -I/home/ycao/opt/src/threadpool/boost
CXXFLAGS+= -I../include -I/usr/include/boost -I/home/khoran/programs/threadpool-0_2_5-src/threadpool/boost -I/home/khoran/programs/libdb/usr/include/

# change this to help me locate boost library files
#LDFLAGS += -L../src  -L ../../atompair -L/home/ycao/opt64/lib/ -L /home/khoran/programs/libdb/usr/lib/
LDFLAGS += -L../src  -L ../../atompair -L /usr/lib/  # -L /home/khoran/programs/libdb/usr/lib/

LDLIBS += -llshkit -ldescriptors

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
DEPS = $(patsubst %.o,%.d,$(OBJS))
PROGS = $(OBJS:.o=)

.PHONY:	all clean lib

all:	lib $(PROGS)

lib:
	$(MAKE) -C ../src

.cpp.o:
	$(CXX) -c -MMD -MP $(CXXFLAGS) $< -o $@
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

clean:
	rm -f $(OBJS)
	rm -f $(DEPS)
	rm -f $(PROGS)

-include $(DEPS)


