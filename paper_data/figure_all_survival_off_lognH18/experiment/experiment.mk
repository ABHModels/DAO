ALL_P1_DIR = results/all_survival_off_5b7fe321
ALL_P1_SOURCES = $(filter-out maindaocl.cpp source/params.cpp source/cloudy_interface_v2.cpp,$(SRCS)) $(ALL_P1_DIR)/source/maindaocl.cpp $(ALL_P1_DIR)/source/params.cpp $(ALL_P1_DIR)/source/cloudy_interface_v2.cpp
$(ALL_P1_DIR)/maindaocl_all_P1: $(ALL_P1_SOURCES) $(HDRS) $(ALL_P1_DIR)/source/verify_no_survival.h
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.cpp,$^) $(LDFLAGS)
