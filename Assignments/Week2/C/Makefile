HOME       := /u/field/for_robert
BLAS_LIB   := $(HOME)/blis/lib/libblis.a
CBLAS_INC  := $(HOME)/blis/include/blis
FLAME_LIB  := $(HOME)/libflame/lib/libflame.a
FLAME_INC  := $(HOME)/libflame/include

# indicate where the object files are to be created
CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -O3 -Wall -I$(FLAME_INC) -I$(CBLAS_INC) -m64 -msse3 -std=c99 
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.
NFIRST     := 500#     smallest size to be timed.
NLAST      := 6000#    largest size to be timed.
NINC       := 500#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

# indicate where the FLAME include files reside

TEST_OBJS  := driver.o SymMatVec1.o SymMatVec4.o SymMatVec5.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)  $(LAPACK_OBJS)
	$(LINKER) $(TEST_OBJS) $(LAPACK_OBJS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LDFLAGS) -o driver.x

test:   driver.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver.x > output.m

clean:
	rm -f *.o *~ core *.x output.m
