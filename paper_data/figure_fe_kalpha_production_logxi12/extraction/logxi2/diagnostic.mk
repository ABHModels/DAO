FE_LOCAL_DIR = results/fe_lowion_emissivity_effac594_20261009
FE_LOCAL_CPP = $(FE_LOCAL_DIR)/diagnostic.cpp source/cloudy_depth.cpp source/rt_grids.cpp source/radiation.cpp source/params.cpp source/corona_models.cpp source/compton_cross_section.cpp
$(FE_LOCAL_DIR)/diagnostic: $(FE_LOCAL_CPP) source/cloudy_interface_v2.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(FE_LOCAL_CPP) $(LDFLAGS)
