#ifndef PELDOR_CALCULATOR_H
#define PELDOR_CALCULATOR_H

#include "definitions.h"
#include "spin.h"
#include <vector>

class PeldorCalculator
{
private:
	int n_fields;

public:
	PeldorCalculator(int const& num_avg);

	// Simulate the random orientation of a magnetic field
	void initialize_field();

	// Calculate the excitation probabilities for the spin A
	void calculate_spinA_excitation(std::vector<exp_signal> const& exp_signals,
		                            Spin const& spinA);

	// Calculate the probablity of excitation of a spin by the detection pulses
	double det_excitation(std::vector<double> const& freqs,
		                  exp_signal const& signal_exp,
						  Spin const& spin) const;

	// Calculate the probablity of excitation of a spin by the pump pulse
	double pump_excitation(std::vector<double> const& freqs,
		                   exp_signal const& signal_exp,
						   Spin const& spin) const;

	// Save spectrum of the spin system
	void save_spectrum(std::vector<exp_signal> const& exp_signals,
		               Spin const& spinA, Spin const& spinB,
					   output_parameters const& output_param) const;

	// Calculate the PELDOR signal for a chromosome
	std::vector<double> calculate_peldor_signal(std::vector<double> const& model_param,
		                                        exp_signal const& signal_exp,
												Spin const& spinA, Spin const& spinB,
												std::vector<int> const& param_numbers,
												std::vector<int> const& param_modes,
												std::mt19937& urng,
												int const& n_offset) const;

	// Set the number of a geometric parameter
	double set_param_value(std::vector<double> const& model_param,
		                   int const& param_number_mean,
						   int const& param_number_width,
						   int const& param_mode_width,
						   std::mt19937& urng) const;

	// RMSD
	double rmsd(std::vector<double> signal_calc, exp_signal const& signal_exp) const;

	// RMSD over Pearson Correlation Coefficient
	double rmsd_over_pearson(std::vector<double> signal_calc, exp_signal const& signal_exp) const;

	// Pearson Correlation Coefficient
	double pearson(std::vector<double> signal_calc, exp_signal const& signal_exp) const;

	// Calculate the form-factor of the PELDOR signal for a chromosome
	std::vector<double> calculate_form_factor(std::vector<double> const& model_param,
		                                      exp_signal const& signal_exp,
											  Spin const& spinA, Spin const& spinB,
											  std::vector<int> const& param_numbers,
											  std::vector<int> const& param_modes,
											  std::mt19937& urng, int const& n_offset) const;

	double threshold;
	std::vector<std::vector<double>> field_orient;
	std::vector<double> gFactorsA;
	std::vector<std::vector<double>> detProbA_allSignals;
	std::vector<std::vector<double>> pumpProbA_allSignals;
};

#endif