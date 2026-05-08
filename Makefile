# ============================================================
# Makefile for DAOv2.0 — Compton RT + Cloudy
# ============================================================
#
# Edit the paths below to match your system.
#
# CLOUDY_SRC : path to your built Cloudy source/ directory
#              (must contain libcloudy.a and cloudyconfig.h)
# XSPEC_LIB  : path to HEASoft Xspec libraries, typically
#              $HEADAS/lib after sourcing headas-init.sh
#
# Remember to `source $HEADAS/headas-init.sh` before running
# any production target — the Xspec shared libraries need it.
# ============================================================

# Cloudy paths (adjust for your system)
CLOUDY_SRC  ?= /path/to/cloudy/source
CLOUDY_LIB  ?= $(CLOUDY_SRC)

# Xspec / HEASoft library path (adjust for your system)
XSPEC_LIB   ?= $(HEADAS)/lib

# Compiler
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -Wall \
           -DSYS_CONFIG=\"$(CLOUDY_SRC)/cloudyconfig.h\" \
           -I$(CLOUDY_SRC) \
           -Isource
LDFLAGS  = -L$(CLOUDY_LIB) -lcloudy -L$(XSPEC_LIB) -lXSFunctions -lm

# Source files by module
SRCS = maindaocl.cpp \
       source/rt_grids.cpp \
       source/params.cpp \
       source/save_results.cpp \
       source/radiation.cpp \
       source/corona_models.cpp \
       source/cloudy_interface_v2.cpp \
       source/compton_cross_section.cpp \
       source/compton_kernel.cpp \
       source/compton_rt.cpp \
       source/source.cpp \
       source/production.cpp \
       source/test_rt.cpp

OBJS = $(SRCS:.cpp=.o)

# Headers (for dependency tracking)
HDRS = source/rt_grids.h \
       source/params.h source/save_results.h source/constants.h \
       source/radiation.h source/corona_models.h \
       source/cloudy_interface.h source/cloudy_exception.h \
       source/compton_cross_section.h source/compton_kernel.h \
       source/compton_rt.h source/source.h \
       source/production.h source/test_rt.h

# Targets
TARGET = maindaocl

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Standalone kernel tests (no Cloudy dependency)
TEST_KERNEL_FLAGS = -std=c++17 -O3 -Wall -Isource

test_kernel_compare: source/test_kernel_compare.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm

test_kernel_norm: source/test_kernel_norm.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm

test_kernel_sym: source/test_kernel_sym.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm

test_kernel_compare_old: source/test_kernel_compare_old.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm

clean:
	rm -f $(OBJS) $(TARGET) \
	      test_kernel_compare test_kernel_norm test_kernel_sym test_kernel_compare_old \
	      normalize_kernel
