#Define executable name
EXECNAME = astroproject

# Define the standard compiler options
CC	= gcc
OPT	= -fomit-frame-pointer -O2
DEBUG	= -g 
PROF	= -pg
LDFLAGS	= -lm -lgslcblas -Wall -lgsl

# Define sourcefiles and objectfiles
SRC	= *.c *.h
OBJ	= *.o

# Complete compilation of the program, optimized
$(EXECNAME): *.c *.h
	$(CC) $(OPT) $(SRC) -o $(EXECNAME) $(LDFLAGS) 
	/bin/rm -rf *.o

# Debug compilation, no optimization or parallelization by default
debug:
	$(CC) $(DEBUG) $(SRC) $(LDFLAGS) -o $(EXECNAME) 

prof:
	$(CC) $(PROF) $(SRC) $(LDFLAGS) -o $(EXECNAME) 

# Compilation of the documentation
doc:	$(SRC)
	doxygen comment_project
