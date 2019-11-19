#For GNU compiler
CXX=g++
CXXFLAGS=-Wall -Wconversion -O3

#For PGI compiler
#CXX=pgc++
#CXXFLAGS=-O3 --gnu

#For CAPS compiler
#CXX=hmpp g++
#CXXFLAGS=-Wall -Wconversion -O3

RM=rm -f
LDFLAGS=-g
LDLIBS=-lfftw3 -lfftw3f -lm -lrt

DESTS=experiment generate_graphs kaiserbessel
SRCS:=$(filter-out $(DESTS:%=%.cc), $(wildcard *.cc))
OBJS:=$(subst .cc,.o,$(SRCS))

all: $(DESTS)

$(DESTS): $(OBJS)

clean:
	$(RM) *.o $(DESTS)
