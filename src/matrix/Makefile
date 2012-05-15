.phony: clean all
EXEC = format read
CXXFLAGS = -I ../atompair
LDFLAGS = -static 
OBJECTS = ../atompair/profiling.o
.cc.o:
	${CXX} -c ${CXXFLAGS} $<

all: ${EXEC}

format: format.o ${OBJECTS}
	${CXX} -o $@ $@.o ${OBJECTS} ${LDFLAGS}

read: read.o ${OBJECTS}
	${CXX} -o $@ $@.o ${OBJECTS} ${LDFLAGS}

clean:
	rm -f *.o ${EXEC}
