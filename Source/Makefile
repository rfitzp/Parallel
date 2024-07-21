# Makefile for Program PARALLEL

F_COMPILER = gfortran
F_FLAGS    =
F_INCLUDE  = -I/opt/local/include

C_COMPILER = clang++ -O3 -std=c++17
C_FLAGS    = -w 
C_INCLUDE  = -I/opt/local/include 
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas -lnetcdf

F_SOURCE = NameListRead.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

C_HEADER = Parallel.h Faddeeva.h
C_SOURCE = main.cpp Parallel.cpp Faddeeva.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../parallel

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(F_OBJECT) $(EXECUTABLE)

