#include "spin.h"
#include <cmath>

Spin::Spin() {}

void Spin::initialize()
{
	if (n.empty()) {
		Nstates = 1;
		Ncomp = 1;
		Icomp.reserve(Ncomp);
		Icomp.push_back(1);
	}
	else if (n.size() == 1) {
		Nstates = static_cast<size_t>(pow((2 * I[0] + 1), n[0]));
		Ncomp = static_cast<size_t>(2 * I[0] * n[0] + 1);
		Icomp.reserve(Ncomp);
		int I_index = static_cast<int>(2 * I[0] - 1);
		int n_index = n[0] - 1;
		int intensity;
		for (size_t i = 0; i < Ncomp; i++) {
			intensity = esr_weights[I_index][n_index][i];
			Icomp.push_back(intensity);
		} 
	}
	else if (n.size() == 2) {
		Nstates = static_cast<size_t>( pow((2 * I[0] + 1), n[0]) * pow((2 * I[1] + 1), n[1]) );
		size_t Ncomp1 = static_cast<size_t>(2 * I[0] * n[0] + 1);
		size_t Ncomp2 = static_cast<size_t>(2 *I[1] * n[1] + 1 );		
		Ncomp = Ncomp1 * Ncomp2;
		Icomp.reserve(Ncomp);
		int I_index1 = static_cast<int>(2 * I[0] - 1);
		int n_index1 = n[0] - 1;
		int I_index2 = static_cast<int>(2 * I[1] - 1);
		int n_index2 = n[1] - 1;
		int intensity;
		for (size_t i = 0; i < Ncomp1; i++) {
			for (size_t j = 0; j < Ncomp2; j++) {
				intensity = esr_weights[I_index1][n_index1][i] * esr_weights[I_index2][n_index2][j];
				Icomp.push_back(intensity);
			}
		}
	}
}

double Spin::calculate_gfactor(std::vector<double> const& proj) const 
{
	double gfactor_eff = sqrt(g[0]*g[0]*proj[0]*proj[0] + 
		                      g[1]*g[1]*proj[1]*proj[1] + 
							  g[2]*g[2]*proj[2]*proj[2]);
	return gfactor_eff;
}

std::vector<double> Spin::calculate_resfreq(std::vector<double> const& proj, 
                                            exp_signal const& signal_exp,
                                            std::mt19937& urng) const
{
	std::vector<double> resfreq; resfreq.reserve(Ncomp);
	double freq(0);
	// A normal distribution
	std::normal_distribution<double> normal_distr(0.0, 1.0);
	// Calculate g-factor
	double gfactor_eff(0), dgfactor_eff(0);
	gfactor_eff = sqrt(g[0]*g[0]*proj[0]*proj[0] + 
		               g[1]*g[1]*proj[1]*proj[1] + 
					   g[2]*g[2]*proj[2]*proj[2]);
	if (gStrain) {
		dgfactor_eff = (g[0]*gStrain[0]*proj[0]*proj[0] + 
			            g[1]*gStrain[1]*proj[1]*proj[1] +
						g[2]*gStrain[2]*proj[2]*proj[2]) / gfactor_eff;
		gfactor_eff += fwhm2sd * dgfactor_eff * normal_distr(urng);
	}
	// Inhomogenious linewidth
	double dfreq(0);
	// Calculate hyperfine coupling
	// The case of 0 nuclei
	if (n.empty()) {
		double dfreq = 1e-3 * 0.5 * lwpp * normal_distr(urng);
		freq = F * gfactor_eff * signal_exp.magnField + dfreq;
		resfreq.push_back(freq);
	}
	// The case of 1 sort of nuclei
	if (n.size() == 1) {
		double A_eff(0), dA_eff(0);
		A_eff = 1e-3 * sqrt(A[0]*A[0]*proj[0]*proj[0] +
			                A[1]*A[1]*proj[1]*proj[1] +
							A[2]*A[2]*proj[2]*proj[2]);
		if (AStrain) {
			dA_eff = 1e-6 * (A[0]*AStrain[0]*proj[0]*proj[0] +
				             A[1]*AStrain[1]*proj[1]*proj[1] +
							 A[2]*AStrain[2]*proj[2]*proj[2]) / A_eff;
			A_eff += fwhm2sd * dA_eff * normal_distr(urng);
		}
		double m = -I[0] * n[0]; // Nuclear quantum number
		for (size_t i = 0; i < Ncomp; ++i) {
			dfreq = 1e-3 * 0.5 * lwpp * normal_distr(urng);
			freq = F * gfactor_eff * signal_exp.magnField + A_eff * m + dfreq;
			resfreq.push_back(freq);
			++m;
		}
	}
	// The case of 2 sorts of nuclei
	if (n.size() == 2) {
		// Hyperfine coupling constants (GHz)
		double A_eff1(0), A_eff2(0), dA_eff(0);
		A_eff1 = 1e-3 * sqrt(A[0]*A[0]*proj[0]*proj[0] + 
			                 A[1]*A[1]*proj[1]*proj[1] + 
							 A[2]*A[2]*proj[2]*proj[2]);
		A_eff2 = 1e-3 * sqrt(A[3]*A[3]*proj[0]*proj[0] + 
			                 A[4]*A[4]*proj[1]*proj[1] + 
							 A[5]*A[5]*proj[2]*proj[2]);
		if (AStrain) {
			dA_eff = 1e-6 * (A[0]*AStrain[0]*proj[0]*proj[0] + 
				             A[1]*AStrain[1]*proj[1]*proj[1] + 
							 A[2]*AStrain[2]*proj[2]*proj[2]) / A_eff1;
			A_eff1 += fwhm2sd * dA_eff * normal_distr(urng);
		}
		// Number of spectral components for each nucleus
		size_t Ncomp1 = static_cast<size_t>(2 * I[0] * n[0] + 1);
		size_t Ncomp2 = static_cast<size_t>(2 * I[1] * n[1] + 1);
		// Nuclear quantum numbers
		double m1 = -I[0] * n[0];
		double m2 = -I[1] * n[1];
		for (size_t i = 0; i < Ncomp1; ++i) {
			m2 = -I[1] * n[1];
			for (size_t j = 0; j < Ncomp2; ++j) {
				dfreq = 1e-3 * 0.5 * lwpp * normal_distr(urng);
				freq = F * gfactor_eff * signal_exp.magnField + A_eff1 * m1  + A_eff2 * m2 + dfreq;
				resfreq.push_back(freq);
				++m2;
			}
			++m1;
		}
	}
	return resfreq;
}