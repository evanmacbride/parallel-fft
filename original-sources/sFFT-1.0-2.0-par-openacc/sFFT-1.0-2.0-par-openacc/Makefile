#For GNU compiler
ifeq ($(CC),gcc)
	CC=gcc 
	CXX=g++ 
	CFLAGS=  -O3 -g -fopenmp -std=gnu99 -Wall -ffast-math -march=native -I/glb/home/uscwc6/include
	CXXFLAGS= -O3 -g -fopenmp -Wall -ffast-math -march=native -I/glb/home/uscwc6/include
	LDFLAGS= -O3 -g -fopenmp -ffast-math -march=native -L/glb/home/uscwc6/lib
	LDLIBS=-lfftw3 -lm -lfftw3_threads -lpthread -lstdc++ #-lfftw3f 
endif

#For Intel compiler
ifeq ($(CC),icc)
	CC=icc
	CXX=icpc
	CFLAGS=  -std=gnu99 -O3 -openmp -fp-model fast=1 -finline -ipo -I/glb/home/uscwc6/include#-no-prec-div
	CXXFLAGS= -O3  -fp-model fast=1 -finline -ipo  #-no-prec-div
	LDFLAGS=  -O3 -openmp  -fp-model fast=1 -finline -ipo -L/glb/home/uscwc6/lib#-no-prec-div -finline
	LDLIBS= -lm -lfftw3 -lfftw3_threads
endif

#For MA
ifeq ($(OPT),ma)
	CFLAGS= -std=gnu99 -O3 -openmp  -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDLIBS= -O3 -openmp -L/glb/data/comp_mod/Tools/fftw/lib -lfftw3 -lfftw3_threads -lm 
endif

#For MA
ifeq ($(OPT),ma)
	CFLAGS= -std=gnu99 -O3 -openmp -g -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDLIBS= -O3 -openmp -L/glb/data/comp_mod/Tools/fftw/lib -lfftw3 -lfftw3_threads -lm 
	CXXFLAGS= -O3  -fp-model fast=1 -finline -openmp -g -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDFLAGS= -O3 -openmp  -fp-model fast=1 -finline #-ipo -no-prec-div -finline
endif

#For PGI compiler
ifeq ($(CC),pgcc)
	PGI_CUDA_PATH =/opt/pgi/linux86-64/2013/cuda/4.2
	FFTW_PATH =/glb/home/uscwc6
	CC=pgcc
	CXX=pgCC
	CFLAGS= -I$(FFTW_PATH)/include -I$(PGI_CUDA_PATH)/include -c99 -O3 -mp -fast -Mipa=fast,inline -Minline  -acc -Minfo=accel -ta=nvidia:fastmath,keepgpu,time,cc35 -Mcuda
	CXXFLAGS= -c99 -O3 -fast -Mipa=fast,inline -Minline
	LDFLAGS= -mp -acc
	LDLIBS= -L$(FFTW_PATH)/lib -L$(PGI_CUDA_PATH)/lib64 -lfftw3 -lm -lrt -lfftw3_threads -lpthread -lstdc++ -mp -acc -lcudart -lcufft 
endif

#For HMPP compiler
ifeq ($(OPT), hmpp)
  CC=hmpp gcc
  CXX = g++
  CFLAGS=  --nvcc-options -Xptxas=-v,-arch=sm_20,-fmad=false -O3 -g -fopenmp -std=gnu99 -Wall -ffast-math -march=native -I/glb/home/uscwc6/include
  CXXFLAGS= -O3 -g -fopenmp -Wall -ffast-math -march=native
  LDFLAGS= -O3 -g -fopenmp -ffast-math -march=native
  LDLIBS=-lfftw3 -lm -lfftw3_threads -lpthread -lstdc++

endif

#RM=rm -f

# SRCS=experiment.c computefourier.c fftw.c filters.c timer.c utils.c wrapper.cc
# DESTS=experiment.o computefourier.o fftw.o filters.o timer.o wrapper.o utils.o

C_SRCS       = $(wildcard *.c)
C_OBJS       = $(C_SRCS:.c=.o)
CXX_SRCS     = $(wildcard *.cc)
CXX_OBJS     = $(CXX_SRCS:.cc=.o)

SRCS         = $(C_SRCS) $(CXX_SRCS)
OBJS         = $(C_OBJS) $(CXX_OBJS)

.PHONY: all	clean

experiment: $(OBJS)
	$(CC) $^ -o $@ $(LDLIBS) 

all:	experiment

clean:
	$(RM) $(OBJS) *.gpu ./experiment
