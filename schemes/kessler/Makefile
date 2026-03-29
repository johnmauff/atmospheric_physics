# ============================================================
# Minimal Fortran Makefile (explicit ordering, no magic)
# ============================================================

# ընտր compiler: gfortran or nvhpc
#COMPILER ?= gfortran
COMPILER ?= nvhpc

# ------------------------------------------------------------
# Compiler selection
# ------------------------------------------------------------

ifeq ($(COMPILER),gnu)
  FC = gfortran
  FFLAGS = -O3 -Wall -cpp
endif

ifeq ($(COMPILER),nvhpc)
  FC = nvfortran
  # FFLAGS = -O3 -Mpreprocess -acc=host -Mnofma
  # FFLAGS = -O3 -Mpreprocess -acc -mp=gpu -gpu=mem:separate,math_uniform -DUSE_GPU
   FFLAGS = -O3 -Mpreprocess -acc -mp=gpu -gpu=mem:separate,math_uniform,nofma -Minfo=acc -DUSE_GPU
endif

# ------------------------------------------------------------
# Target
# ------------------------------------------------------------

TARGET = test_kessler

# ------------------------------------------------------------
# Build (EXPLICIT ORDER — THIS IS THE KEY)
# ------------------------------------------------------------

all: $(TARGET)

compare: 
	./test_kessler >& test.$(COMPILER).log
	diff test.$(COMPILER).log golden-nofma.$(COMPILER).log

$(TARGET): ccpp_kinds.o kessler.o driver_kessler.o
	$(FC) $(FFLAGS) -o $(TARGET) \
	    ccpp_kinds.o kessler.o driver_kessler.o

# Compile in strict dependency order

ccpp_kinds.o: ccpp_kinds.F90
	$(FC) $(FFLAGS) -c ccpp_kinds.F90

kessler.o: kessler.F90 ccpp_kinds.o
	$(FC) $(FFLAGS) -c kessler.F90

driver_kessler.o: driver_kessler.F90 kessler.o
	$(FC) $(FFLAGS) -c driver_kessler.F90

# ------------------------------------------------------------
# تنظيف
# ------------------------------------------------------------

clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: all clean
