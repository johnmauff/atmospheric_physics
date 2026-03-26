# ============================================================
# GNU Fortran Makefile for Kessler test driver
# ============================================================

# Compiler
FC       = gfortran

# Executable name
TARGET   = test_kessler

# Source files
SRCS     = kessler.F90 \
           test_kessler_driver.F90

# Object files (automatic)
OBJS     = $(SRCS:.F90=.o)

# Default build type
BUILD    ?= release

# ------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------

ifeq ($(BUILD),debug)
  FFLAGS = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS = -O3 -march=native -Wall
endif

# Preprocessing enabled because files use .F90
CPPFLAGS = -cpp

# ------------------------------------------------------------
# Build rules
# ------------------------------------------------------------

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

# ------------------------------------------------------------
# Utility targets
# ------------------------------------------------------------

debug:
	$(MAKE) BUILD=debug

release:
	$(MAKE) BUILD=release

clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: all clean debug release
