# Makefile for Program PARALLEL

C_COMPILER = clang++
C_FLAGS    = -w -O3 -std=c++17
C_INCLUDE  = -I/opt/local/include 
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgsl -lgslcblas -lnetcdf

C_HEADER = Parallel.h Faddeeva.h
C_SOURCE = main.cpp Parallel.cpp Faddeeva.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../parallel

$(EXECUTABLE):	 	$(C_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

