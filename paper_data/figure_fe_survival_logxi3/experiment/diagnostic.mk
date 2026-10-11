FE_CONTROL_DIR = results/fe_k_survival_a518ab75/transfer_control_20261009
$(FE_CONTROL_DIR)/transfer_control: $(FE_CONTROL_DIR)/transfer_control.cpp source/cloudy_depth.o source/rt_grids.o source/radiation.o source/params.o source/corona_models.o source/compton_cross_section.o source/compton_rt.o source/thermal_balance.o source/avg_compton_kernel.o source/compton_kernel.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
