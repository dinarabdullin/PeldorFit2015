#include "read_config.h"
#ifdef WIN32
#pragma comment(lib, "libconfig")
#endif
#include "libconfig.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

void read_config(const char* path_config, 
	             std::vector<exp_signal>& signals_exp,
				 Spin& spinA, Spin& spinB, 
				 std::vector<int>& param_numbers,
				 std::vector<int>& param_modes,
				 std::vector<bound>& param_bounds,
				 genetic_parameters& genetic_param,
				 output_parameters& output_param,
				 errorplot_parameters& errorplot_param)
{
	// Read a configuration file
	config_t cfg, *cf;
	cf = &cfg;
	config_init(cf);
	if (!config_read_file(cf, path_config)) {
		std::cerr << config_error_file(cf) << ": " << config_error_line(cf) << " - " << config_error_text(cf) << std::endl;
		config_destroy(cf);
		exit(EXIT_FAILURE);
	}
	// Enable the auto convertion of data types
	config_set_auto_convert(cf, true);
	////////////////////////////////////////////////////////////////////////////
	//                        EXPERIMENTAL PARAMETERS                         //
	////////////////////////////////////////////////////////////////////////////
	const config_setting_t *experimentals = config_lookup(cf, "experimentals");
	size_t numOfSignals = config_setting_length(experimentals);
	std::cout << "  Number of PELDOR signals is " << numOfSignals << std::endl;
	const config_setting_t *set;
	const char *path = NULL;
	std::string dataline;
	double timeval, signalval;
	std::vector<double> timearray, signalarray;
	double detPiLength_temp, detPiHalfLength_temp, pumpPiLength_temp, detFreq_temp, pumpFreq_temp, magnField_temp;
	exp_signal signal_temp;
	signals_exp.reserve(numOfSignals);
	for (size_t i = 0; i < numOfSignals; ++i) {
		set = config_setting_get_elem(experimentals, i);
		path = config_setting_get_string(config_setting_get_member(set, "filename"));
		std::ifstream datafile(path, std::ios::in);
		if (!datafile.is_open()) {
			std::cerr << " Error: Could not open the data file!" << std::endl;
			exit(EXIT_FAILURE);
		}
		else {
			while (getline(datafile, dataline)) {
				std::stringstream datastream(dataline);
				datastream >> timeval >> signalval;
				timearray.push_back(timeval);
				signalarray.push_back(signalval);
			}
			signal_temp.timeAxis = timearray;
			signal_temp.signalAxis = signalarray;
			timearray.clear();
			signalarray.clear();
		}
		detPiLength_temp = config_setting_get_float(config_setting_get_member(set, "detPiLength"));
		detPiHalfLength_temp = config_setting_get_float(config_setting_get_member(set, "detPiHalfLength"));
		pumpPiLength_temp = config_setting_get_float(config_setting_get_member(set, "pumpPiLength"));
		detFreq_temp = config_setting_get_float(config_setting_get_member(set, "detFreq"));
		pumpFreq_temp = config_setting_get_float(config_setting_get_member(set, "pumpFreq"));
		magnField_temp = config_setting_get_float(config_setting_get_member(set, "magnField"));
		signal_temp.detPiHalfLength = detPiHalfLength_temp;
		signal_temp.detPiLength = detPiLength_temp;
		signal_temp.pumpPiLength = pumpPiLength_temp;
		signal_temp.detFreq = detFreq_temp;
		signal_temp.pumpFreq = pumpFreq_temp;
		signal_temp.magnField = magnField_temp;
		signals_exp.push_back(signal_temp);
	}
	////////////////////////////////////////////////////////////////////////////
	//								 SPIN PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	// Read the spectroscopic parameters of spin A
	const config_setting_t *pspinA = config_lookup(cf, "spinA");
	for (size_t i = 0; i < 3; ++i) {
		spinA.g[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "g"), i);
	}
	size_t sizeA = config_setting_length(config_setting_get_member(pspinA, "n"));
	for (size_t i = 0; i < sizeA; ++i) {
		spinA.n.push_back(config_setting_get_int_elem(config_setting_get_member(pspinA, "n"), i));
		spinA.I.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "I"), i));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 0));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 1));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 2));
	}
	if (config_setting_length(config_setting_get_member(pspinA, "gStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i) 
			spinA.gStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "gStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinA.gStrain[i] = 0.0;
	}
	if (config_setting_length(config_setting_get_member(pspinA, "AStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i) 
			spinA.AStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "AStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinA.AStrain[i] = 0.0;
	}
	spinA.lwpp = config_setting_get_float(config_setting_get_member(pspinA, "lwpp"));
	spinA.initialize();
	// Read the spectroscopic parameters of spin B
	const config_setting_t *pspinB = config_lookup(cf, "spinB");
	for (size_t i = 0; i < 3; ++i) {
		spinB.g[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "g"), i);
	}
	size_t sizeB = config_setting_length(config_setting_get_member(pspinB, "n"));
	for (size_t i = 0; i < sizeB; ++i) {
		spinB.n.push_back(config_setting_get_int_elem(config_setting_get_member(pspinB, "n"), i));
		spinB.I.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "I"), i));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 0));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 1));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 2));
	}
	if (config_setting_length(config_setting_get_member(pspinB, "gStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i)
			spinB.gStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "gStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinB.gStrain[i] = 0.0;
	}
	if (config_setting_length(config_setting_get_member(pspinB, "AStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i)
			spinB.AStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "AStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinB.AStrain[i] = 0.0;
	}
	spinB.lwpp = config_setting_get_float(config_setting_get_member(pspinB, "lwpp"));
	spinB.initialize();
	////////////////////////////////////////////////////////////////////////////
	//							  FITTING PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	const config_setting_t *fit_param = config_lookup(cf, "parameters");
	const config_setting_t *single_param;
	int fit_enable(0);
	int mode(0);
	bound single_bound;
	int count(0);
	for (int i = 0; i < 28; ++i) {
		// Read parameter after parameter
		single_param = config_setting_get_elem(fit_param, i);
		// Check if this parameter is enabled
		fit_enable = config_setting_get_int(config_setting_get_member(single_param, "opt"));
		if (fit_enable) {
			// Associate a model parameter with the index of the corresponding fitting parameter
			param_numbers.push_back(count);
			++count;
			// Read out the distribution types of the fitting parameters
			mode = config_setting_get_int(config_setting_get_member(single_param, "mode"));
			param_modes.push_back(mode);
			// Read out the bounds for the fitting parameters
			single_bound.lower = config_setting_get_float_elem(config_setting_get_member(single_param, "range"), 0);
			single_bound.upper = config_setting_get_float_elem(config_setting_get_member(single_param, "range"), 1);
			if (std::count(angular_indices, angular_indices + 20, i)) {
				single_bound.lower *= deg2rad;
				single_bound.upper *= deg2rad;
			}
			param_bounds.push_back(single_bound);
			if ((i == 27) && (mode == 1)) {
				for (size_t j = 0; j < (numOfSignals - 1); ++j) {
					param_bounds.push_back(single_bound);
					++count;
				}
			}
		}
		else {
			param_numbers.push_back(-1);
			param_modes.push_back(0);
		}
	}
	param_bounds.reserve(count);
	////////////////////////////////////////////////////////////////////////////
	//			   	   PARAMETERS OF THE GENETIC ALGORITHM                    //
	////////////////////////////////////////////////////////////////////////////
	config_setting_t *genetic = config_lookup(cf, "genetic");
	genetic_param.num_generations_max = config_setting_get_int(config_setting_get_member(genetic, "num_generations_max"));
	genetic_param.size_generation = config_setting_get_int(config_setting_get_member(genetic, "size_generation"));
	genetic_param.prob_crossover = config_setting_get_float(config_setting_get_member(genetic, "prob_crossover"));
    genetic_param.prob_mutation = config_setting_get_float(config_setting_get_member(genetic, "prob_mutation"));
	genetic_param.merit_function = config_setting_get_int(config_setting_get_member(genetic, "merit_function"));
	genetic_param.num_avg = config_setting_get_int(config_setting_get_member(genetic, "num_averages"));
	////////////////////////////////////////////////////////////////////////////
	//				               OUTPUT PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	config_setting_t *output = config_lookup(cf, "output");
	const char *pfileroot = NULL;
	pfileroot = config_setting_get_string(config_setting_get_member(output, "directory"));
	std::string output_dir(pfileroot);
	output_param.directory = output_dir;
    output_param.record_spectrum = config_setting_get_int(config_setting_get_member(output, "record_spectrum"));
    output_param.record_score = config_setting_get_int(config_setting_get_member(output, "record_score"));
    output_param.record_parameters = config_setting_get_int(config_setting_get_member(output, "record_parameters"));
    output_param.record_fit = config_setting_get_int(config_setting_get_member(output, "record_fit"));
    output_param.record_form_factor = config_setting_get_int(config_setting_get_member(output, "record_form_factor"));
    output_param.record_symmetric_solutions = config_setting_get_int(config_setting_get_member(output, "record_symmetric_solutions"));
    output_param.record_error_plot = config_setting_get_int(config_setting_get_member(output, "record_error_plot"));
    config_setting_t *ep_var = config_setting_get_member(output, "error_plot_variables");
	int nErrorPlots = config_setting_length(ep_var);
	int nVarParam(0);
	if (nErrorPlots == 0) {
		output_param.record_error_plot = 0;
	}
	else {
		output_param.error_plot_variables.reserve(nErrorPlots);
		for (int i = 0; i < nErrorPlots; ++i) {
			nVarParam = config_setting_length(config_setting_get_elem(ep_var, i));
			std::vector<int> varParam;  varParam.reserve(nVarParam);
			for (int j = 0; j < nVarParam; ++j) varParam.push_back(config_setting_get_int_elem(config_setting_get_elem(ep_var,i),j));
			output_param.error_plot_variables.push_back(varParam);
			varParam.clear();
		}
	}
	output_param.error_plot_size = config_setting_get_int(config_setting_get_member(output, "error_plot_size"));
	////////////////////////////////////////////////////////////////////////////
	//				            'ERROR PLOT ONLY' MODE                        //
	////////////////////////////////////////////////////////////////////////////
	const config_setting_t *errorplot = config_lookup(cf, "error_plot_only");
	errorplot_param.enable = config_setting_get_int(config_setting_get_member(errorplot, "enable"));
	if (errorplot_param.enable) {
		config_setting_t *ep_var = config_setting_get_member(errorplot, "error_plot_variables");
		int nErrorPlots = config_setting_length(ep_var);
		int nVarParam(0);
		if (nErrorPlots != 0) {
			errorplot_param.error_plot_variables.reserve(nErrorPlots);
			for (int i = 0; i < nErrorPlots; ++i) {
				nVarParam = config_setting_length(config_setting_get_elem(ep_var, i));
				std::vector<int> varParam;  varParam.reserve(nVarParam);
				for (int j = 0; j < nVarParam; ++j) {
					varParam.push_back(config_setting_get_int_elem(config_setting_get_elem(ep_var,i),j));
				}
				errorplot_param.error_plot_variables.push_back(varParam);
				varParam.clear();
			}
		}
		errorplot_param.error_plot_size = config_setting_get_int(config_setting_get_member(errorplot, "error_plot_size"));
		// Read out the file with the optimized parameters
		errorplot_param.optimized_genes.reserve(param_bounds.size());
		// Set the number of lines to read
		int n_datalines(28);
		if ((param_numbers[27] != -1) && (param_modes[27] == 1)) n_datalines += (numOfSignals - 1);
		// Set the current number of line (-1 line corresponds to the capitel)
		int count(-1); 
		// Read the file line by line
		double datavalue;
		std::string dataline;
		const char *pInputFile = NULL;
		pInputFile = config_setting_get_string(config_setting_get_member(errorplot, "input_directory"));
		std::ifstream inputFile(pInputFile, std::ios::in);
		if (!inputFile.is_open()) {
			std::cerr << " Error: Could not open the data file!" << std::endl;
			exit(EXIT_FAILURE);
		}
		else {
			while (getline(inputFile, dataline)) {
				if ((count >= 0) && (count < n_datalines) && (param_numbers[count] != -1)) {
					std::string subline = dataline.substr(30);
					std::stringstream datastream(subline);
					datastream >> datavalue;
					if (std::count(angular_indices, angular_indices + 20, count)) {
						errorplot_param.optimized_genes.push_back(datavalue * deg2rad);
					}
					else {
						errorplot_param.optimized_genes.push_back(datavalue);
					}
				}
				++count;
			}
		}
		const char *pOutputFile = NULL;
		pOutputFile = config_setting_get_string(config_setting_get_member(errorplot, "output_directory"));
		std::string outputFile(pOutputFile);
		errorplot_param.output_directory = outputFile;
	}
	// Close a configuration file
	config_destroy(cf);
};