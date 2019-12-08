#ifndef SPIN_H
#define SPIN_H

#include "definitions.h"
#include <vector>
#include <random>

class Spin
{
public:
	double g[3];
	std::vector<int> n;
	std::vector<double> I;
	std::vector<double> A;
	double gStrain[3];
	double AStrain[3];
	double lwpp;
	size_t Nstates; // Number of microstates
	size_t Ncomp; // Number of spectral components
	std::vector<int> Icomp; // Intensities of the individual spectral components

	Spin();

	// Calculate number and intensities of ESR lines
	void initialize();

	// Compute g-factor
	double calculate_gfactor(std::vector<double> const& proj) const;

	// Compute resonance frequencies
	std::vector<double> calculate_resfreq(std::vector<double> const& proj, 
		                                  exp_signal const& signal_exp,
										  std::mt19937& urng) const;
};

#endif