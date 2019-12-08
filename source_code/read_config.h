#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include "definitions.h"
#include "spin.h"
#include <vector>

void read_config(const char* path_config, 
	             std::vector<exp_signal>& signals_exp,
				 Spin& spinA, Spin& spinB, 
				 std::vector<int>& param_numbers,
				 std::vector<int>& param_modes,
				 std::vector<bound>& param_bounds,
				 genetic_parameters& genetic_param,
				 output_parameters& output_param,
				 errorplot_parameters& errorplot_param);
#endif