OBJDIR = obj
SRCDIR = src
BINDIR = bin
LIBDIR = src/lib

SRC = $(wildcard $(SRCDIR)/*.c)
OBJ = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRC))
BIN = $(filter-out $(FILTER), $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%, $(SRC)))
LIB = $(wildcard $(LIBDIR)/*.c)
##############################################################################
OS = $(shell uname)

ifeq ($(OS),Darwin)
$(info Detected MacOS as Operating System. Omitting OpenMP compile flags!)
CC = clang
CFLAGS = -std=c99 -O2 -march=native -fopenmp
endif

ifeq ($(OS),Linux)
$(info Detected Linux as Operating System.)
CC = gcc
CFLAGS = -std=c99 -O2 -march=native -fopenmp
endif

MASTER_EQ_DEPS 	= $(OBJDIR)/masterEquation.o $(OBJDIR)/graphlib.o $(OBJDIR)/matgen.o  $(OBJDIR)/matvec.o $(OBJDIR)/solv.o $(OBJDIR)/io.o
MF_FULL_DEPS	= $(OBJDIR)/mf_full.o $(OBJDIR)/io.o $(OBJDIR)/matvec.o
SS_GMP_DEPS		= $(OBJDIR)/steadystates_gmp.o $(OBJDIR)/io.o

.PHONY = debug warn asan no-gmp no-omp verbose all dirs
.DEFAULT = all

ifdef debug
CFLAGS += -ggdb
CFLAGS:=$(filter-out -O2,$(CFLAGS))
endif

ifdef warn
CFLAGS += -Wall -Wextra -Wpedantic -Wformat -Wstrict-aliasing -Wmain
CFLAGS += -Wpointer-arith -Wunused-parameter -Wno-uninitialized
CFLAGS += -pedantic -Wmisleading-indentation -Wmissing-include-dirs
endif

ifdef asan
CFLAGS += -fsanitize=address,undefined -static-libasan
LDFLAGS += -fsanitize=address,undefined
endif
CFLAGS += -D_DEFAULT_SOURCE -MMD -MP

ifeq ($(OS), Linux)
LDFLAGS += -lm -ldl -lgmp
endif

ifdef no-gmp
LDFLAGS:=$(filter-out -lgmp,$(LDFLAGS))
SRC:=$(filter-out $(SRCDIR)/steadystates_gmp.c,$(SRC))
BIN:=$(filter-out $(BINDIR)/steadystates_gmp,  $(BIN))
OBJ:=$(filter-out $(OBJDIR)/steadystates_gmp.o,$(OBJ))
endif

ifdef no-omp
CFLAGS:=$(filter-out -fopenmp,$(CFLAGS))
endif

ifdef verbose
CFLAGS += -v
LDFLAGS += -v
endif

all: dirs $(OBJ) $(BIN)

dirs:
	mkdir -p ./bin ./obj


#Linking
$(BINDIR)/masterEquation: $(MASTER_EQ_DEPS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

ifdef no-gmp
$(BINDIR)/steadystates_gmp: $(SS_GMP_DEPS)
	$(info GMP support disabled)
else
$(BINDIR)/steadystates_gmp: $(SS_GMP_DEPS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@
	
endif

$(BINDIR)/mf_full: $(MF_FULL_DEPS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@
	
#Assembly
$(OBJDIR)/%.o: $(LIBDIR)/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c 
	$(CC) $(CFLAGS) $^ -c -o $@ $(LDFLAGS)



clean:
	rm -f bin/* obj/*.o log/make.log
