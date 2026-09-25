# ============================================================
# Makefile for Compton RT + Cloudy
# ============================================================

# Cloudy paths — the ONLY thing you change per version
CLOUDY_SRC  = /path/to/cloudy/source
#eg. CLOUDY_SRC=  ./cloudyc25/source
CLOUDY_LIB  = $(CLOUDY_SRC)
CLOUDY_ROOT = $(CLOUDY_SRC)/..

# Compiler
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -Wall -pthread \
           -DSYS_CONFIG=\"$(CLOUDY_SRC)/cloudyconfig.h\" \
           -I$(CLOUDY_SRC) \
           -Isource

XSPEC_LIB = /Users/ym.huang/Downloads/heasoft-6.33.2/aarch64-apple-darwin22.3.0/lib

# Auto-detect VectorHash: present in C25+, absent in older Cloudy.
# Picks lib64 if it exists, else lib32, else nothing.
VH_LIB := $(firstword $(wildcard $(CLOUDY_ROOT)/library/vectorhash/lib64 \
                                 $(CLOUDY_ROOT)/library/vectorhash/lib32))
ifneq ($(VH_LIB),)
  VH_FLAGS = -L$(VH_LIB) -lvhsum -Wl,-rpath,$(VH_LIB)
endif

# -lvhsum MUST come after -lcloudy (libcloudy references VectorHash)
LDFLAGS  = -L$(CLOUDY_LIB) -lcloudy $(VH_FLAGS) -L$(XSPEC_LIB) -lXSFunctions -lm

# Source files by module
SRCS = maindaocl.cpp \
       source/rt_grids.cpp \
       source/params.cpp \
       source/save_results.cpp \
       source/radiation.cpp \
       source/corona_models.cpp \
       source/cloudy_interface_v2.cpp \
       source/cloudy_depth.cpp \
       source/compton_cross_section.cpp \
       source/compton_kernel.cpp \
       source/avg_compton_kernel.cpp \
       source/compton_rt.cpp \
       source/source.cpp \
       source/production.cpp \
       source/test_rt.cpp

OBJS = $(SRCS:.cpp=.o)

# Headers (for dependency tracking)
HDRS = source/rt_grids.h \
       source/params.h source/save_results.h source/constants.h \
       source/radiation.h source/corona_models.h \
       source/cloudy_interface.h source/cloudy_exception.h source/cloudy_depth.h \
       source/compton_cross_section.h source/compton_kernel.h source/avg_compton_kernel.h \
       source/kernel_payload.h source/kernel_row_spool.h source/rt_parallel.h \
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

# Standalone kernel test (no Cloudy dependency)
TEST_KERNEL_FLAGS = -std=c++17 -O3 -Wall -pthread -Isource

test_kernel_storage: source/test_kernel_storage.cpp source/kernel_payload.h source/kernel_row_spool.h source/rt_parallel.h
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $<

test_kernel_norm: source/test_kernel_norm.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm
multiscat: source/multiscat.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm
dump_kernel_slice: source/dump_kernel_slice.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm
clean:
	rm -f $(OBJS) $(TARGET) test_kernel_norm test_kernel_storage multiscat dump_kernel_slice
