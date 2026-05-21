MF = Makefile

CC = g++

CFLAGS = -D_USE_64 -msse3 -O3 -g -fomit-frame-pointer -funroll-loops
LFLAGS = -std=c++17 -O3 -DNDEBUG -I ./libsdsl/include/ -L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 -Wl,-rpath=$(PWD)/libsdsl/lib

EXE = lamas
SRC = lamas.cc
OBJ = $(SRC:.cc=.o)

HD = lamas.h Makefile

.SUFFIXES:
.SUFFIXES: .cc .o

# Compile .cc -> .o
.cc.o:
	$(CC) $(CFLAGS) $(LFLAGS) -c $< -o $@

# Default target
all: $(EXE)

# Link executable
$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LFLAGS)

# Dependencies
$(OBJ): $(MF) $(HD)

# Clean object files and executable
clean:
	rm -f $(OBJ) $(EXE) *~

clean-all:
	rm -f $(OBJ) $(EXE) *~
