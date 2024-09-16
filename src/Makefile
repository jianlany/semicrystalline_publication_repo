# C++ compiler (and linker)
CXX=ccache g++

# Include libraries
STRICT_WARNINGS=-Wall -Wno-unused-function -Wno-sign-compare #-Werror 
CFLAGS+= -I. -DDEBUGGING_ON ${STRICT_WARNINGS}
LFLAGS+=

# Extracts the current commit (short) hash.
gitrev=${shell git log --pretty=format:'%h' -n 1}
# Lists all modified files (should be none for production runs.)
is_modified=${shell git diff-index --name-only HEAD .}
CFLAGS+=-DGIT_REVISION=\"${gitrev}\"
CFLAGS+=-DMODIFIED_FILES="\"${is_modified}\""

LFLAGS+=-llapack -lpthread -lrt 
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=molbuilder

all: $(SRC) $(EXE) $(OBJ) tags
all: CFLAGS+=-O3 -fopenmp -march=native
all: LFLAGS+= -fopenmp

debug: $(SRC) $(EXE) $(OBJ)
debug: CFLAGS+=-g -fopenmp
debug: LFLAGS+=-g -fopenmp

profile: $(SRC) $(EXE) $(OBJ)
profile: CFLAGS+=-O2 -pg
profile: LFLAGS+=-pg

doc:
	@cd ../doc && doxygen

tags: $(SRC)
	@ctags $(wildcard *.h */*.h *.cpp */*.cpp)

$(EXE): $(OBJ) $(EX_LIB_FILE)
	$(CXX) $(OBJ) $(LFLAGS) -o $@

.cpp.o:
	$(CXX) -c -std=c++0x $(CFLAGS) $< -o $@

clean:
	@rm -f $(OBJ) tags

main.o: .FORCE
.FORCE:
