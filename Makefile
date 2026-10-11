# ============================================================
# Makefile for Compton RT + Cloudy
# ============================================================

# Cloudy paths — the ONLY thing you change per version
CLOUDY_SRC  = /Users/xe26734/Downloads/cloudy/source
#eg. CLOUDY_SRC=  ./cloudyc25/source
CLOUDY_LIB  = $(CLOUDY_SRC)
CLOUDY_ROOT = $(CLOUDY_SRC)/..

# Compiler
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -Wall -pthread \
           -DSYS_CONFIG=\"$(CLOUDY_SRC)/cloudyconfig.h\" \
           -I$(CLOUDY_SRC) \
           -Isource

XSPEC_LIB = /Users/xe26734/Downloads/heasoft-6.37.1/aarch64-apple-darwin25.6.0/lib

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
       source/thermal_balance.cpp \
       source/thermal_production.cpp \
       source/test_rt.cpp

OBJS = $(SRCS:.cpp=.o)

# Headers (for dependency tracking)
HDRS = source/rt_grids.h \
       source/params.h source/save_results.h source/constants.h source/run_log.h \
       source/radiation.h source/corona_models.h \
       source/cloudy_interface.h source/cloudy_exception.h source/cloudy_depth.h \
       source/compton_cross_section.h source/compton_kernel.h source/avg_compton_kernel.h \
       source/kernel_payload.h source/kernel_row_spool.h source/rt_parallel.h \
       source/kernel_quadrature.h \
       source/compton_rt.h source/cell_transfer.h source/incidence_boundary.h source/source.h \
       source/production.h source/test_rt.h source/thermal_balance.h source/thermal_production.h

# Targets
TARGET = maindaocl

.PHONY: all clean

all: $(TARGET)

# Reproduce compPS's five outgoing quadrature nodes without changing the
# production executable or its default eight-angle grid. Same sources/solver.
.PHONY: compps_native_benchmark
compps_native_benchmark: results/compps_native/maindaocl

results/compps_native/maindaocl: $(SRCS) $(HDRS)
	mkdir -p results/compps_native
	$(CXX) $(CXXFLAGS) -DDAO_RT_ANGLES=10 -o $@ $(SRCS) $(LDFLAGS)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Standalone kernel test (no Cloudy dependency)
TEST_KERNEL_FLAGS = -std=c++17 -O3 -Wall -pthread -Isource

test_compton_transfer: source/test_compton_transfer.o source/compton_rt.o source/thermal_balance.o source/compton_kernel.o source/avg_compton_kernel.o source/compton_cross_section.o source/rt_grids.o source/radiation.o source/corona_models.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_kernel_storage: source/test_kernel_storage.cpp source/kernel_payload.h source/kernel_row_spool.h source/rt_parallel.h
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $<

# Local-only tests are available when their ignored sources are present.
ifneq ($(wildcard source/test_incidence_boundary.cpp),)
test_incidence_boundary: source/test_incidence_boundary.cpp source/params.cpp source/rt_grids.cpp source/radiation.cpp source/corona_models.cpp source/incidence_boundary.h
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.cpp,$^) $(LDFLAGS)
endif

ifneq ($(wildcard source/test_cloudy_depth.cpp),)
test_cloudy_depth: source/test_cloudy_depth.o source/save_results.o source/cloudy_depth.o source/cloudy_interface_v2.o source/compton_cross_section.o source/rt_grids.o source/radiation.o source/corona_models.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
endif

ifneq ($(wildcard source/test_ion_fractions.cpp),)
test_ion_fractions: source/test_ion_fractions.o source/save_results.o source/cloudy_depth.o source/cloudy_interface_v2.o source/compton_cross_section.o source/rt_grids.o source/radiation.o source/corona_models.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
endif

test_scattering_density: source/test_scattering_density.cpp source/source.cpp source/rt_grids.cpp source/compton_kernel.cpp source/avg_compton_kernel.cpp source/compton_cross_section.cpp $(HDRS)
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $(filter %.cpp,$^) -lm

test_kernel_low_temp: source/test_kernel_low_temp.cpp source/source.cpp source/rt_grids.cpp source/compton_kernel.cpp source/avg_compton_kernel.cpp source/compton_cross_section.cpp $(HDRS)
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $(filter %.cpp,$^) -lm
multiscat: source/multiscat.cpp source/rt_grids.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm
dump_kernel_slice: source/dump_kernel_slice.cpp source/compton_kernel.cpp source/compton_cross_section.cpp
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $^ -lm
clean:
	rm -f test_ion_fractions source/test_ion_fractions.o
	rm -f $(OBJS) $(TARGET) test_compton_transfer source/test_compton_transfer.o test_incidence_boundary test_kernel_low_temp test_kernel_storage test_scattering_density test_cloudy_depth source/test_cloudy_depth.o test_thermal_balance test_thermal_directional test_thermal_cloudy source/test_thermal_cloudy.o multiscat dump_kernel_slice

# Thermal closure regression: conservation, detailed balance, root failures.
test_thermal_balance: source/test_thermal_balance.cpp source/thermal_balance.cpp source/avg_compton_kernel.cpp source/compton_kernel.cpp source/compton_cross_section.cpp source/rt_grids.cpp source/test_local_response.h $(HDRS)
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $(filter %.cpp,$^) -lm

ifneq ($(wildcard source/test_thermal_directional.cpp),)
test_thermal_directional: source/test_thermal_directional.cpp source/thermal_balance.cpp source/avg_compton_kernel.cpp source/compton_kernel.cpp source/compton_cross_section.cpp source/rt_grids.cpp source/test_local_response.h $(HDRS)
	$(CXX) $(TEST_KERNEL_FLAGS) -o $@ $(filter %.cpp,$^) -lm
endif

ifneq ($(wildcard source/test_thermal_cloudy.cpp),)
test_thermal_cloudy: source/test_thermal_cloudy.o source/compton_rt.o source/production.o source/thermal_production.o source/thermal_balance.o source/save_results.o source/cloudy_depth.o source/cloudy_interface_v2.o source/compton_cross_section.o source/avg_compton_kernel.o source/compton_kernel.o source/rt_grids.o source/radiation.o source/corona_models.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
endif
