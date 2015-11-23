#include "peldor_calculator.h"
#include "rotations.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <functional>
#include <algorithm>
#include <array>

PeldorCalculator::PeldorCalculator(int const& num_avg) 
{
	threshold = 0.001;
	n_fields = num_avg;
	initialize_field();
}

void PeldorCalculator::initialize_field()
{
	field_orient.reserve(n_fields);
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 engine(seq);
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	// Calculate random orientations of the magnetic field vector
	double fxi(0), fphi(0);
	std::vector<double> single_orient;
	single_orient.reserve(3);
	for (int f = 0; f < n_fields; ++f) {
		fphi = 0.0 + (2.0*PI - 0.0) * uniform_distr(engine);
		fxi = acos(uniform_distr(engine));
		if (uniform_distr(engine) < 0.5) fxi = PI - fxi;
		single_orient.push_back(sin(fxi)*cos(fphi));
		single_orient.push_back(sin(fxi)*sin(fphi));
		single_orient.push_back(cos(fxi));
		field_orient.push_back(single_orient);
		single_orient.clear();
	}
}

void PeldorCalculator::calculate_spinA_excitation(std::vector<exp_signal> const& signals_exp,
	                                              Spin const& spinA)
{
	gFactorsA.reserve(n_fields);
	size_t numOfSignals = signals_exp.size();
	detProbA_allSignals.reserve(numOfSignals);
	pumpProbA_allSignals.reserve(numOfSignals);
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 engine(seq);
	// Calculate the excitation probabilities for the spin A
	double gfactorA(0);
	std::vector<double> resfreqA; resfreqA.reserve(spinA.Ncomp);
	double detProbA(0), pumpProbA(0);
	std::vector<double> detProbA_oneSignal; detProbA_oneSignal.reserve(n_fields);
	std::vector<double> pumpProbA_oneSignal; pumpProbA_oneSignal.reserve(n_fields);
	for (size_t i = 0; i < numOfSignals; ++i) {
		for (int f = 0; f < n_fields; ++f) {
			// Compute an effective g-factor of the spin A
			gfactorA = spinA.calculate_gfactor(field_orient[f]);
			// Compute resonance frequencies of the spin A
			resfreqA = spinA.calculate_resfreq(field_orient[f], signals_exp[i], engine);
			// Compute the probability of spin A to be excited by detection pulses
			detProbA = det_excitation(resfreqA, signals_exp[i], spinA);
			// Compute the probability of spin A to be excited by a pump pulse
			pumpProbA = pump_excitation(resfreqA, signals_exp[i], spinA);
			// Store the calculated values
			gFactorsA.push_back(gfactorA);
			detProbA_oneSignal.push_back(detProbA);
			pumpProbA_oneSignal.push_back(pumpProbA);
			resfreqA.clear();
		}
		detProbA_allSignals.push_back(detProbA_oneSignal);
		pumpProbA_allSignals.push_back(pumpProbA_oneSignal);
		detProbA_oneSignal.clear();
		pumpProbA_oneSignal.clear();
	}
}

double PeldorCalculator::det_excitation(std::vector<double> const& freqs,
	                                    exp_signal const& signal_exp,
										Spin const& spin) const
{
	double detProb(0);
	double detPiHalfBW = 0.25 / signal_exp.detPiHalfLength;
	double detPiBW = 0.50 / signal_exp.detPiLength;
	double freqEff(0), freqEff2(0), prob(0);
	for (size_t k = 0; k < spin.Ncomp; k++) {
		if (detPiHalfBW == detPiBW) {
			freqEff = sqrt(pow(signal_exp.detFreq - freqs[k], 2) + pow(detPiHalfBW, 2));
			prob = pow((detPiHalfBW / freqEff) * sin(2*PI * freqEff * signal_exp.detPiHalfLength), 5);
		}
		else {
			freqEff = sqrt(pow(signal_exp.detFreq - freqs[k], 2) + pow(detPiHalfBW, 2));
			freqEff2 = sqrt(pow(signal_exp.detFreq - freqs[k], 2) + pow(detPiBW, 2));
			prob = (detPiHalfBW / freqEff) * sin(2*PI * freqEff * signal_exp.detPiHalfLength) * 
				   pow((detPiBW / freqEff2) * sin(0.5 * 2*PI * freqEff2 * signal_exp.detPiLength), 4);
		}
		detProb += fabs(static_cast<double>(spin.Icomp[k]) * prob);
	}
	detProb /= static_cast<double>(spin.Nstates);
	return detProb;
}

double PeldorCalculator::pump_excitation(std::vector<double> const& freqs,
	                                     exp_signal const& signal_exp,
										 Spin const& spin) const
{
	double pumpProb(0);
	double pumpPiBW = 0.5 / signal_exp.pumpPiLength;
	double freqEff(0), prob(0);
	for (size_t k = 0; k < spin.Ncomp; k++) {
		freqEff = sqrt(pow(signal_exp.pumpFreq - freqs[k], 2) + pow(pumpPiBW, 2));
		prob = pow((pumpPiBW / freqEff) * sin(0.5 * 2*PI * freqEff * signal_exp.pumpPiLength), 2);
		pumpProb += fabs(static_cast<double>(spin.Icomp[k]) * prob);
	}
	pumpProb /= static_cast<double>(spin.Nstates);
	return pumpProb;
}

void PeldorCalculator::save_spectrum(std::vector<exp_signal> const& signals_exp,
	                                 Spin const& spinA, Spin const& spinB,
									 output_parameters const& output_param) const
{
	// Set the file in which the the spectrum will be recorded
	std::ostringstream filename;
	filename << output_param.directory << "spectrum.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 engine(seq);
	// Determine the lowest and the highest frequencies
	double freqMin, freqMax, freqMinA, freqMaxA, freqMinB, freqMaxB;
	std::vector<double> resfreqA; resfreqA.reserve(spinA.Ncomp);
	std::vector<double> resfreqB; resfreqB.reserve(spinB.Ncomp);
	for (int f = 0; f < n_fields; ++f) {
		// Compute the resonance frequences of the spin A
		resfreqA = spinA.calculate_resfreq(field_orient[f], signals_exp[0], engine);
		// Compute the resonance frequence of the spin B
		resfreqB = spinB.calculate_resfreq(field_orient[f], signals_exp[0], engine);
		// Search for the lowest and the highest frequencies
		freqMinA = *std::min_element(resfreqA.begin(),resfreqA.end());
		freqMaxA = *std::max_element(resfreqA.begin(),resfreqA.end());
		freqMinB = *std::min_element(resfreqB.begin(),resfreqB.end());
		freqMaxB = *std::max_element(resfreqB.begin(),resfreqB.end());
		if (f == 0) {
			if (freqMinA < freqMinB) freqMin = freqMinA;
			else freqMin = freqMinB;
			if (freqMaxA > freqMaxB) freqMax = freqMaxA;
			else freqMax = freqMaxB;
		}
		else {
			if (freqMinA < freqMin) freqMin = freqMinA;
			if (freqMinB < freqMin) freqMin = freqMinB;
			if (freqMaxA > freqMax) freqMax = freqMaxA;
			if (freqMaxB > freqMax) freqMax = freqMaxB;
		}
		resfreqA.clear();
		resfreqB.clear();
	}
	// Increase slightly the frequency range
	freqMin -= 0.100;
	freqMax += 0.100;
	// Create the axes of the spectrum
	double df = 0.001; // step of 1 MHz
	freqMin = ceil(freqMin/df) * df; // round to the 1 MHz precision
	freqMax = ceil(freqMax/df) * df; // round to the 1 MHz precision
	int nf = static_cast<int>((freqMax - freqMin) / df); // number of intervals 
	std::vector<double> freqAxis; freqAxis.reserve(nf);  // x-axis
	std::vector<double> specAxis; specAxis.reserve(nf); // y-axis
	int index;
	int index_min = static_cast<int>(freqMin/df) + 1; // Correlate frequencies with the indeces
	// Initialize axes
	double freq(0);
	for (int i = 0; i < nf; i++) {
		freq = freqMin + 0.5 * df * (2*i + 1);
		freqAxis.push_back(freq); // initialize x-axis
		specAxis.push_back(0); // initialize y-axis
	}
	// Calculate the spectrum
	for (int f = 0; f < n_fields; ++f) {
		// Compute the resonance frequences of the spin A
		resfreqA = spinA.calculate_resfreq(field_orient[f], signals_exp[0], engine);
		// Compute the resonance frequence of the spin B
		resfreqB = spinB.calculate_resfreq(field_orient[f], signals_exp[0], engine);
		// Add the calculated frequencies to the spectrum
		for (size_t k = 0; k < spinA.Ncomp; k++) {
			index = static_cast<int>(resfreqA[k]/df) - index_min;
			specAxis[index] += static_cast<double>(spinA.Icomp[k]) / static_cast<double>(spinA.Nstates);
		}
		for (size_t k = 0; k < spinB.Ncomp; k++) {
			index = static_cast<int>(resfreqB[k]/df) - index_min;
			specAxis[index] += static_cast<double>(spinB.Icomp[k]) / static_cast<double>(spinB.Nstates);
		}
		resfreqA.clear();
		resfreqB.clear();
	}
	// Normilize the spectrum and write it into the file
	double specAxisMax = *max_element(specAxis.begin(), specAxis.end());
	file << std::left << std::setw(20) << "Frequency(GHz)" << std::setw(20) << "Amplitude" << std::endl;
	for (int i = 0; i < nf; i++) {
		file << std::left << std::setprecision(5)
			 << std::setw(20) << freqAxis[i] 
		     << std::setw(20) << specAxis[i]/specAxisMax << std::endl;
	}
	file.close();	
}

std::vector<double> PeldorCalculator::calculate_peldor_signal(std::vector<double> const& model_param,
	                                                          exp_signal const& signal_exp,
															  Spin const& spinA, Spin const& spinB,
															  std::vector<int> const& param_numbers,
															  std::vector<int> const& param_modes,
															  std::mt19937& urng,
															  int const& n_offset) const
{
	size_t const n_points = signal_exp.timeAxis.size();
	double* signal_calc_array = new double[n_points];
	for (size_t t = 0; t < n_points; ++t) signal_calc_array[t] = 0.0;
	// Ititialize some variables
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	double dist(0), xi(0), phi(0), alpha(0), betta(0), gamma(0);
	double jc(0);
	std::vector<double> projB; projB.reserve(3);
	double gFactorB(0);
	std::vector<double> resfreqB; resfreqB.reserve(spinB.Ncomp);
	double detProbA(0), pumpProbA(0), detProbB(0), pumpProbB(0);
	double signal_amplitude(0);
	double cosDipAngle(0);
	double freq_dd(0), freq_dd_constant(0), freq_dd_angular(0);
	// Set the probability of the first mean distance with respect to the second mean distance
	double ratio(0);
	if (param_numbers[24] == -1) ratio = 1.0;
	else                         ratio = model_param[param_numbers[24]];
	// Set the iversion efficiency of the pump pulse
	double pumpEfficiency(1);
	if (param_numbers[27] == -1) pumpEfficiency = 1;
	else {
		if (param_modes[27] == 0) pumpEfficiency = model_param[param_numbers[27]];
		else if (param_modes[27] == 1) pumpEfficiency = model_param[param_numbers[27]+n_offset];
	}
	// Iterate over orientations of the field and relative orientations of the spins
	for (int f = 0; f < n_fields; ++f) {
		// Set the values of the geometric parameters and the j coupling
		if (uniform_distr(urng) <= ratio) {
			dist  = set_param_value(model_param,param_numbers[0],param_numbers[1],param_modes[1],urng);
			xi    = set_param_value(model_param,param_numbers[2],param_numbers[3],param_modes[3],urng);
			phi   = set_param_value(model_param,param_numbers[4],param_numbers[5],param_modes[5],urng);
			alpha = set_param_value(model_param,param_numbers[6],param_numbers[7],param_modes[7],urng);
			betta = set_param_value(model_param,param_numbers[8],param_numbers[9],param_modes[9],urng);
			gamma = set_param_value(model_param,param_numbers[10],param_numbers[11],param_modes[11],urng);
			jc    = set_param_value(model_param,param_numbers[25],param_numbers[26],param_modes[26],urng);
		}
		else {
			dist  = set_param_value(model_param,param_numbers[12],param_numbers[13],param_modes[1],urng);
			xi    = set_param_value(model_param,param_numbers[14],param_numbers[15],param_modes[3],urng);
			phi   = set_param_value(model_param,param_numbers[16],param_numbers[17],param_modes[5],urng);
			alpha = set_param_value(model_param,param_numbers[18],param_numbers[19],param_modes[7],urng);
			betta = set_param_value(model_param,param_numbers[20],param_numbers[21],param_modes[9],urng);
			gamma = set_param_value(model_param,param_numbers[22],param_numbers[23],param_modes[11],urng);
			jc    = set_param_value(model_param,param_numbers[25],param_numbers[26],param_modes[26],urng);
		}
		// Rotation matrix between the spin A and spin B coordinate systems
		ZXZEulerRM const* rotationMatrix = new ZXZEulerRM(alpha,betta,gamma);
		// Compute the orientation of the magnetic field in the frame of the spin B
		projB = rotationMatrix->rotate_vector(field_orient[f]);
		// Compute the effective g-factor of the spin B
		gFactorB = spinB.calculate_gfactor(projB);
		// Compute resonance frequencies of the spin B
		resfreqB = spinB.calculate_resfreq(projB,signal_exp,urng);
		// Compute the probability of spin B to be excited by detection pulses
		detProbB = det_excitation(resfreqB,signal_exp,spinB);
		// Read out the excitation probabilities for the spin A
		detProbA = detProbA_allSignals[n_offset][f];
		pumpProbA = pumpProbA_allSignals[n_offset][f];
		// Calculate the amplitude of the PELDOR signal
		signal_amplitude += (detProbA + detProbB);
		if (detProbA > threshold) {
			// Compute the probability of spin B to be excited by a pump pulse
			pumpProbB = pump_excitation(resfreqB,signal_exp,spinB);
			// Check if the excitation of the spins by the pump pulses is big enough
			if ((pumpProbA > threshold) || (pumpProbB > threshold)) {
				// Compute the angular part of the dipolar frequency
				cosDipAngle = field_orient[f][0] * cos(phi) * sin(xi) + 
					          field_orient[f][1] * sin(phi) * sin(xi) + 
							  field_orient[f][2] * cos(xi);
				freq_dd_angular = 1.0 - 3.0 * cosDipAngle * cosDipAngle;
				// Compute the dipolar constant
				freq_dd_constant = 52.04 * gFactorsA[f] * gFactorB /(2.0023 * 2.0023 * pow(dist,3));
				// Compute the dipolar frequency
				freq_dd = 2*PI * (freq_dd_constant * freq_dd_angular + jc);
				// Compute the oscillating part of the PELDOR signal
				for (size_t t = 0; t < n_points; ++t) {
					signal_calc_array[t] += pumpEfficiency * (detProbA * pumpProbB + detProbB * pumpProbA) * (1.0 - cos(freq_dd * signal_exp.timeAxis[t]));
				}
			}
		}
		else if ((detProbA < threshold) && (detProbB > threshold) && (pumpProbA > threshold)) {
			// Compute the angular part of the dipolar frequency
			cosDipAngle = field_orient[f][0] * cos(phi) * sin(xi) +
				          field_orient[f][1] * sin(phi) * sin(xi) +
						  field_orient[f][2] * cos(xi);
			freq_dd_angular = 1.0 - 3.0 * cosDipAngle * cosDipAngle;
			// Compute the dipolar constant
			freq_dd_constant = 52.04 * gFactorsA[f] * gFactorB /(2.0023 * 2.0023 * pow(dist,3));
			// Compute the dipolar frequency
			freq_dd = 2*PI * (freq_dd_constant * freq_dd_angular + jc);
			// Compute the oscillating part of the PELDOR signal
			for (size_t t = 0; t < n_points; ++t) {
				signal_calc_array[t] += pumpEfficiency * (detProbB * pumpProbA) * (1.0 - cos(freq_dd * signal_exp.timeAxis[t]));
			}
		}
		// Clean up
		delete rotationMatrix;
		projB.clear();
		resfreqB.clear();
	}
	// Calculate a PELDOR signal, normalize it and save as a vector
	std::vector<double> signal_calc; signal_calc.reserve(n_points);
	double norm = 1.0 / signal_amplitude;
	for (size_t t = 0; t < n_points; ++t) {
		signal_calc_array[t] = (signal_amplitude - signal_calc_array[t]) * norm;
		signal_calc.push_back(signal_calc_array[t]);
	}
	delete [] signal_calc_array; signal_calc_array = nullptr;
	return signal_calc;
}

double PeldorCalculator::set_param_value(std::vector<double> const& model_param,
	                                     int const& param_number_mean,
										 int const& param_number_width,
										 int const& param_mode_width,
										 std::mt19937& urng) const
{
	double param_value(0), param_mean(0), param_width(0);
	if (param_number_mean == -1) {
		param_value = 0;
	}
	else {
		param_mean = model_param[param_number_mean];
		if (param_number_width == -1) {
			param_value = param_mean;
		}
		else {
			param_width = model_param[param_number_width];
			if (param_mode_width == 0) {
				std::uniform_real_distribution<> uniform_distr(param_mean-0.5*param_width, param_mean+0.5*param_width);
				param_value = uniform_distr(urng);
			}
			else if (param_mode_width == 1) {
				std::normal_distribution<> normal_distr(param_mean, param_width);
				param_value = normal_distr(urng);
			}
		}
	}
	return param_value;
}

double PeldorCalculator::rmsd(std::vector<double> signal_calc, exp_signal const& signal_exp) const
{
	double rmsd(0);
	size_t const n_points = signal_exp.signalAxis.size();
	for (size_t t = 0; t < n_points; ++t) rmsd += pow((signal_calc[t] - signal_exp.signalAxis[t]), 2);
	rmsd = sqrt(rmsd / static_cast<double>(n_points));
	return rmsd;
}

double PeldorCalculator::rmsd_over_pearson(std::vector<double> signal_calc, exp_signal const& signal_exp) const
{
	size_t const n_points = signal_exp.signalAxis.size();
	double signal_exp_mean(0), signal_exp_variance(0);
	double signal_calc_mean(0), signal_calc_variance(0);
	double rmsd(0);
	double pcc(0);
	double rmsd_pcc(0);
	for (size_t t = 0; t < n_points; ++t) {
		signal_exp_mean += signal_exp.signalAxis[t];
		signal_calc_mean  += signal_calc[t];
	}
	signal_exp_mean /= static_cast<double>(n_points);
	signal_calc_mean /= static_cast<double>(n_points);
	for (size_t t = 0; t < n_points; ++t) {
		signal_exp_variance += (signal_exp.signalAxis[t] - signal_exp_mean) * (signal_exp.signalAxis[t] - signal_exp_mean);
		signal_calc_variance += (signal_calc[t] - signal_calc_mean) * (signal_calc[t] - signal_calc_mean);
		rmsd += (signal_calc[t] - signal_exp.signalAxis[t]) * (signal_calc[t] - signal_exp.signalAxis[t]);
		pcc += (signal_calc[t] - signal_calc_mean) * (signal_exp.signalAxis[t] - signal_exp_mean);
	}
	pcc /= sqrt(signal_calc_variance * signal_exp_variance);
	rmsd = sqrt(rmsd / static_cast<double>(n_points));
	rmsd_pcc = rmsd / pcc;
	return rmsd_pcc;
}

double PeldorCalculator::pearson(std::vector<double> signal_calc, exp_signal const& signal_exp) const
{
	size_t const n_points = signal_exp.signalAxis.size();
	double signal_exp_mean(0), signal_exp_variance(0);
	double signal_calc_mean(0), signal_calc_variance(0);
	double pcc(0);
	for (size_t t = 0; t < n_points; ++t) {
		signal_exp_mean += signal_exp.signalAxis[t];
		signal_calc_mean  += signal_calc[t];
	}
	signal_exp_mean /= static_cast<double>(n_points);
	signal_calc_mean /= static_cast<double>(n_points);
	for (size_t t = 0; t < n_points; ++t) {
		signal_exp_variance += (signal_exp.signalAxis[t] - signal_exp_mean) * (signal_exp.signalAxis[t] - signal_exp_mean);
		signal_calc_variance += (signal_calc[t] - signal_calc_mean) * (signal_calc[t] - signal_calc_mean);
		pcc += (signal_calc[t] - signal_calc_mean) * (signal_exp.signalAxis[t] - signal_exp_mean);
	}
	pcc /= sqrt(signal_calc_variance * signal_exp_variance);
	return pcc;
}

std::vector<double> PeldorCalculator::calculate_form_factor(std::vector<double> const& model_param,
	                                                        exp_signal const& signal_exp,
															Spin const& spinA, Spin const& spinB,
															std::vector<int> const& param_numbers,
															std::vector<int> const& param_modes,
															std::mt19937& urng, int const& n_offset) const
{
	double form_factor_array[90] = {0};
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	double dist(0), xi(0), phi(0), alpha(0), betta(0), gamma(0);
	std::vector<double> projB; projB.reserve(3);
	double gFactorB(0);
	std::vector<double> resfreqB; resfreqB.reserve(spinB.Ncomp);
	double detProbA(0), pumpProbA(0), detProbB(0), pumpProbB(0);
	double dipAngle(0), cosDipAngle(0);
	int dipAngle_index(0);
	// Set the probability of the first mean distance with respect to the second mean distance
	double ratio(0);
	if (param_numbers[24] == -1) ratio = 1.0;
	else                         ratio = model_param[param_numbers[24]];
	// Iterate over orientations of the field and relative orientations of the spins
	for (int f = 0; f < n_fields; ++f) {
		// Set the values of the geometric parameters and the j coupling
		if (uniform_distr(urng) <= ratio) {
			dist  = set_param_value(model_param,param_numbers[0],param_numbers[1],param_modes[1],urng);
			xi    = set_param_value(model_param,param_numbers[2],param_numbers[3],param_modes[3],urng);
			phi   = set_param_value(model_param,param_numbers[4],param_numbers[5],param_modes[5],urng);
			alpha = set_param_value(model_param,param_numbers[6],param_numbers[7],param_modes[7],urng);
			betta = set_param_value(model_param,param_numbers[8],param_numbers[9],param_modes[9],urng);
			gamma = set_param_value(model_param,param_numbers[10],param_numbers[11],param_modes[11],urng);
		}
		else {
			dist  = set_param_value(model_param,param_numbers[12],param_numbers[13],param_modes[1],urng);
			xi    = set_param_value(model_param,param_numbers[14],param_numbers[15],param_modes[3],urng);
			phi   = set_param_value(model_param,param_numbers[16],param_numbers[17],param_modes[5],urng);
			alpha = set_param_value(model_param,param_numbers[18],param_numbers[19],param_modes[7],urng);
			betta = set_param_value(model_param,param_numbers[20],param_numbers[21],param_modes[9],urng);
			gamma = set_param_value(model_param,param_numbers[22],param_numbers[23],param_modes[11],urng);
		}
		// Rotation matrix between the spin A and spin B coordinate systems
		ZXZEulerRM const* rotationMatrix = new ZXZEulerRM(alpha,betta,gamma);
		// Compute the orientation of the magnetic field in the frame of the spin B
		projB = rotationMatrix->rotate_vector(field_orient[f]);
		// Compute the effective g-factor of the spin B
		gFactorB = spinB.calculate_gfactor(projB);
		// Compute resonance frequencies of the spin B
		resfreqB = spinB.calculate_resfreq(projB,signal_exp,urng);
		// Compute the probability of spin B to be excited by the detection pulses
		detProbB = det_excitation(resfreqB,signal_exp,spinB);
		// Compute the probability of spin B to be excited by the pump pulse
		pumpProbB = pump_excitation(resfreqB,signal_exp,spinB);
		// Read out the excitation probabilities for the spin A
		detProbA = detProbA_allSignals[n_offset][f];
		pumpProbA = pumpProbA_allSignals[n_offset][f];
		// Calculate the dipolar angle
		cosDipAngle = field_orient[f][0] * cos(phi) * sin(xi) +
			          field_orient[f][1] * sin(phi) * sin(xi) +
					  field_orient[f][2] * cos(xi);
		dipAngle = acos(cosDipAngle) * rad2deg;
		if (dipAngle < 0)  dipAngle = -dipAngle;
		if (dipAngle > 90) dipAngle = 180 - dipAngle;
		dipAngle_index = static_cast<int>( std::floor(dipAngle) );
		if (dipAngle_index > 89) dipAngle_index = 89;
		// Relate the calculated probabilities to the corresponding dipolar angle
		form_factor_array[dipAngle_index] += (detProbA * pumpProbB + detProbB * pumpProbA);
		// Clean up
		delete rotationMatrix;
		projB.clear();
		resfreqB.clear();
	}
	// Normalize the form factor
	std::vector<double> form_factor; form_factor.reserve(90);
	double integral(0);
	for (int i = 0; i < 90; ++i) integral += form_factor_array[i];
	integral *= PI / 180;
	for (int i = 0; i < 90; ++i) form_factor.push_back( form_factor_array[i] / integral);
	return form_factor;
}