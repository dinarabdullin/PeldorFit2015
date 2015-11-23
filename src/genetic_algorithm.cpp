#include "genetic_algorithm.h"
#include "rotations.h"
#include "tbb/tbb.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>
#include <random>
#include <functional>
#include <algorithm>
#include <array>

//------------------------------ CHROMOSOME CLASS ------------------------------

Chromosome::Chromosome(std::vector<bound> const& param_bounds,
					   std::mt19937& urng)
{
	// Set the size of a chromosome
	size = param_bounds.size();
	genes.reserve(size);
	// Create initial genes
	double gene(0);
	for (size_t i = 0; i < size; ++i) {
		gene = create_random_gene(param_bounds[i].lower,param_bounds[i].upper,urng);
		genes.push_back(gene);
	}
}

double Chromosome::create_random_gene(double const& gene_min,
	                                  double const& gene_max,
									  std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(gene_min, gene_max);
	return uniform_distr(urng);
}

//------------------------------ GENERATION CLASS ------------------------------

Generation::Generation(genetic_parameters const& genetic_param)
{
	// Set the size of a generation
	size = genetic_param.size_generation;
	chromosomes.reserve(size);
	offspring.reserve(size);
}

void Generation::create_initial_generation(std::vector<bound> const& param_bounds,
	                                       genetic_parameters const& genetic_param)
{
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 engine(seq);
	// Create chromosomes
	for (size_t i = 0; i < size; ++i) {
		Chromosome* chromosome = new Chromosome(param_bounds,engine);
		chromosomes.push_back(*chromosome);
		delete chromosome;
	}
}


void Generation::score_chromosomes(PeldorCalculator const& peldor_calc,
	                               std::vector<exp_signal> const& signals_exp,
								   Spin const& spinA, Spin const& spinB,
								   std::vector<int> const& param_numbers,
								   std::vector<int> const& param_modes,
								   genetic_parameters const& genetic_param)
{
	const size_t numOfSignals = signals_exp.size();
	tbb::parallel_for(size_t(0), size, [&](size_t i) {
		// Initialize random generator
		std::random_device rd;
		std::array<unsigned int, std::mt19937::state_size> seed_data;
		std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
		std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
		std::mt19937 engine(seq);
		double goodnessOfFit(0);
		for (size_t j = 0; j < numOfSignals; ++j) {
			// Calculate the PELDOR signal for each chromosome and for each frequency offset
			std::vector<double> signal_calc; signal_calc.reserve(signals_exp[j].signalAxis.size());
			signal_calc = peldor_calc.calculate_peldor_signal(chromosomes[i].genes,signals_exp[j],spinA,
				                                              spinB,param_numbers,param_modes,engine,j);
			// Calculate the fitness for each chromosome
			if (genetic_param.merit_function == 1) {
				goodnessOfFit += peldor_calc.rmsd(signal_calc,signals_exp[j]);
			}
			else if (genetic_param.merit_function == 2) {
				goodnessOfFit += peldor_calc.rmsd_over_pearson(signal_calc,signals_exp[j]);
			}
			else if (genetic_param.merit_function == 3) {
				goodnessOfFit += peldor_calc.pearson(signal_calc,signals_exp[j]);
			}
		}
		chromosomes[i].fitness = goodnessOfFit;
	});
}

void Generation::sort_chromosomes()
{ 
	std::sort(chromosomes.begin(), chromosomes.end()); 
}

void Generation::produce_offspring(genetic_parameters const& genetic_param, std::vector<bound> const& param_bounds)
{
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 engine(seq);
	// Check if the number of chromosomes is even
	size_t pairs(0);
	if (size % 2 == 0) {
		pairs = size / 2;
	}
	else {
		pairs = (size + 1) / 2;
	}
	// Select the pairs of parents and produce an offspring
	for (size_t i = 0; i < pairs; ++i) {
		// Select parents via tournament selection
		Chromosome parent1 = chromosomes[tournament_selection(engine)];
		Chromosome parent2 = chromosomes[tournament_selection(engine)];
		// Crossover parents
		crossover_chromosomes(parent1, parent2, genetic_param.prob_crossover, engine);
		// Mutate parents
		mutate_chromosome(parent1, genetic_param.prob_mutation, param_bounds, engine);
		mutate_chromosome(parent2, genetic_param.prob_mutation, param_bounds, engine);
		// Save new chromosomes
		offspring.push_back(parent1);
		if (2 * (i + 1) <= size) {
			offspring.push_back(parent2);
		}
	}
	// Save new generation
	chromosomes = offspring;
	offspring.clear();
}

int Generation::tournament_selection(std::mt19937& urng) const
{
	std::uniform_int_distribution<int> uniform_distr(0, size-1);
	int n1 = uniform_distr(urng);
	int n2 = uniform_distr(urng);
	if (chromosomes[n1].fitness < chromosomes[n2].fitness) return n1;
	else                                                   return n2;
}

void Generation::crossover_chromosomes(Chromosome& chromosome1, Chromosome& chromosome2, 
	                                   double const& prob_crossover, 
									   std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	if (uniform_distr(urng) <= prob_crossover) {
		std::vector<double> genes_temp; genes_temp.reserve(chromosome1.size); // store initial genes of the 1st chromosome
		for (size_t i = 0; i < chromosome1.size; ++i) genes_temp.push_back(chromosome1.genes[i]);
		std::uniform_int_distribution<int> uniform_distr2(1, chromosome1.size - 2); 
		int point = uniform_distr2(urng); // choose a random crossover point
		for (int i = 0; i <= point; ++i) {
			chromosome1.genes[i] = chromosome2.genes[i];
			chromosome2.genes[i] = genes_temp[i];
		}
	}
}

void Generation::mutate_chromosome(Chromosome& chromosome, 
	                               double const& prob_mutation,
								   std::vector<bound> const& param_bounds, 
								   std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	for (size_t i = 0; i < chromosome.size; ++i) {
		if (uniform_distr(urng) <= prob_mutation) {
			chromosome.genes[i] = chromosome.create_random_gene(param_bounds[i].lower, param_bounds[i].upper, urng);
		}
	}
}

//------------------------------ GENETIC ALGORITHM CLASS ------------------------------

GeneticAlgorithm::GeneticAlgorithm() 
{}

void GeneticAlgorithm::run_optimization(std::vector<exp_signal> const& signals_exp,
	                                    Spin const& spinA, Spin const& spinB,
										std::vector<int> const& param_numbers,
										std::vector<int> const& param_modes,
										std::vector<bound> const& param_bounds,
										genetic_parameters const& genetic_param,
										output_parameters const& output_param) const
{
	// Initialize the calculator of PELDOR signals
	PeldorCalculator* peldor_calc = new PeldorCalculator(genetic_param.num_avg);
	peldor_calc->calculate_spinA_excitation(signals_exp, spinA);
	// Calculate the spectrum of the spin system
	if (output_param.record_spectrum) {
		std::cout << "  Recording the spectrum... ";
		peldor_calc->save_spectrum(signals_exp, spinA, spinB, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Start optimization
	int n = 1;
	std::cout << "  Optimization step " << n << "/" << genetic_param.num_generations_max << std::endl;
	// Create a first generation
	Generation* generation = new Generation(genetic_param);
	generation->create_initial_generation(param_bounds, genetic_param);
	// Score first generation
	generation->score_chromosomes(*peldor_calc, signals_exp, spinA, spinB, param_numbers, param_modes, genetic_param);
	generation->sort_chromosomes();
	std::vector<double> best_fitness; best_fitness.reserve(genetic_param.num_generations_max);
	best_fitness.push_back(generation->chromosomes[0].fitness);
	// Proceed with next generations
	for (int n = 2; n <= genetic_param.num_generations_max; ++n) { // loop over generations
		std::cout << "  Optimization step " << n << "/" << genetic_param.num_generations_max << std::endl;
		// Create next generation
		generation->produce_offspring(genetic_param, param_bounds);
		// Score the generation
		generation->score_chromosomes(*peldor_calc, signals_exp, spinA, spinB, param_numbers, param_modes, genetic_param);
		generation->sort_chromosomes();
		best_fitness.push_back(generation->chromosomes[0].fitness);
	}
	// Save the fitness of the best chromosome
	if (output_param.record_score) {
		std::cout << "  Recording the fitness vs the number of optimization steps... ";
		record_score(best_fitness, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Save the calculated fit to the PELDOR signals for the best chromosome
	if (output_param.record_fit) {
		std::cout << "  Recording the fit to the PELDOR signals... ";
		record_fit(*peldor_calc, generation->chromosomes[0], signals_exp, spinA, spinB, param_numbers, param_modes, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Save the genes of the best chromosome
	if (output_param.record_parameters) {
		std::cout << "  Recording the best values of fitting parameters... ";
		record_parameters(generation->chromosomes[0], param_numbers, param_modes, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Record the calculated for-factors for the PELDOR signals
	if (output_param.record_symmetric_solutions) {
		std::cout << "  Recording the form-factors of PELDOR signals... ";
		record_form_factor(*peldor_calc, generation->chromosomes[0], signals_exp, spinA, spinB,
			               param_numbers, param_modes, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Record the symmetry-related sets of fitting parameters
	if (output_param.record_symmetric_solutions) {
		std::cout << "  Recording the symmetry-related sets of fitting parameters... ";
		record_symmetric_parameters(*peldor_calc, generation->chromosomes[0], signals_exp, spinA, spinB, 
			                        param_numbers, param_modes, param_bounds, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Record the error plots
	if (output_param.record_error_plot) {
		std::cout << "  Recording the error plot... ";
		record_error_plot(*peldor_calc, generation->chromosomes[0], signals_exp, spinA, spinB, 
			              param_numbers, param_modes, param_bounds, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Clean up
	delete peldor_calc;
	delete generation;
}

void GeneticAlgorithm::record_fit(PeldorCalculator const& peldor_calc,
	                              Chromosome const& chromosome,
		                          std::vector<exp_signal> const& signals_exp,
								  Spin const& spinA, Spin const& spinB,
								  std::vector<int> const& param_numbers,
								  std::vector<int> const& param_modes,
								  output_parameters const& output_param) const
{
	// Calculate the PELDOR signals
	const size_t numOfSignals = signals_exp.size();
	std::vector<std::vector<double>> signals_calc; signals_calc.reserve(numOfSignals);
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 engine(seq);
	for (size_t j = 0; j < numOfSignals; ++j) {
		// Calculate the PELDOR signal for the best chromosome and different frequency offsets
		std::vector<double> signal_calc; signal_calc.reserve(signals_exp[j].signalAxis.size());
		signal_calc = peldor_calc.calculate_peldor_signal(chromosome.genes, signals_exp[j], spinA, spinB, 
				                                          param_numbers, param_modes, engine, j);
		signals_calc.push_back(signal_calc);
	}
	// Create & open a file
	std::ostringstream filename;
	filename << output_param.directory << "fit.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the names of columns
	size_t const col_width = 10;
	std::ostringstream exp_name;
	std::ostringstream fit_name;
	for (size_t i = 0; i < signals_exp.size(); ++i) {
		exp_name << "exp" << i + 1;
		fit_name << "fit" << i + 1;
		file << std::left << std::setw(col_width) << "t(us)"
			              << std::setw(col_width) << exp_name.str() 
						  << std::setw(col_width) << fit_name.str();
		exp_name.str(std::string()); exp_name.clear();
		fit_name.str(std::string()); fit_name.clear();
	}
	file << std::endl;
	// Calculate the maximal number of time points
	size_t nOfPoints_max(0);
	for (size_t i = 0; i < signals_exp.size(); ++i) {
		if (signals_exp[i].timeAxis.size() > nOfPoints_max) nOfPoints_max = signals_exp[i].timeAxis.size();
	}
	// Write the calculated PELDOR signals and their fits
	for (size_t t = 0; t < nOfPoints_max; ++t) {
		for (size_t i = 0; i < signals_exp.size(); ++i) {
			if (t < signals_exp[i].timeAxis.size()) {
				file << std::left << std::fixed << std::setprecision(5)
					 << std::setw(col_width) << signals_exp[i].timeAxis[t]
					 << std::setw(col_width) << signals_exp[i].signalAxis[t]
					 << std::setw(col_width) << signals_calc[i][t];
			}
			else {
				file << std::setw(3*col_width) << " ";
			}
		}
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_score(std::vector<double> const& best_fitness,
	                                genetic_parameters const& genetic_param,
									output_parameters const& output_param) const
{
	// Create & open a file
	std::ostringstream filename;
	filename << output_param.directory << "score.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the names of columns
	size_t const col_width = 20;
	file << std::left << std::setw(col_width) << "Generation";
	if (genetic_param.merit_function == 1)      file << std::left << std::setw(col_width) << "RMSD" << std::endl;
	else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width) << "RMSD/PCC" << std::endl;
	else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width) << "PCC" << std::endl;
	// Write the number of generation and RMSD
	for (size_t n = 0; n < best_fitness.size(); ++n) {
		file << std::left << std::setw(col_width) << (n + 1);
		file << std::left << std::setw(col_width) << std::fixed << std::setprecision(5) << best_fitness[n];
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_parameters(Chromosome const& chromosome,
		                                 std::vector<int> const& param_numbers,
										 std::vector<int> const& param_modes,
										 output_parameters const& output_param) const
{
	// Create & open a file
	std::ostringstream filename;
	filename << output_param.directory << "parameters.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the names of columns
	size_t const col_width = 30;
	file << std::left << std::setw(col_width) << "Parameter";
	file << std::left << std::setw(col_width) << "Value";
	file << std::endl;
	// Fill in the columns with names of parameters and their values
	const char *param_names[] = { "distance mean 1 (nm)","distance width 1 (nm)","xi mean 1 (deg)","xi width 1 (deg)","phi mean 1 (deg)","phi width 1 (deg)",
		                          "alpha mean 1 (deg)","alpha width 1 (deg)","betta mean 1 (deg)","betta width 1 (deg)","gamma mean 1 (deg)","gamma width 1 (deg)",
								  "distance mean 2 (nm)", "distance width 2 (nm)", "xi mean 2 (deg)", "xi width 2 (deg)", "phi mean 2 (deg)", "phi width 2 (deg)",
								  "alpha mean 2 (deg)", "alpha width 2 (deg)", "betta mean 2 (deg)", "betta width 2 (deg)", "gamma mean 2 (deg)", "gamma width 2 (deg)",
								  "ratio 1 to 2", "J mean (MHz)", "J width (MHz)", "inv. efficiency" };
	for (size_t i = 0; i < 28; ++i) {
		file << std::left << std::setw(col_width) << param_names[i];
		if (param_numbers[i] == -1) {
			file << std::setw(col_width) << " ";
			file << std::endl;
		}
		else {
			if (std::count(angular_indices, angular_indices + 20, i)) {
				file << std::left << std::setw(col_width) << std::fixed << std::setprecision(0) << chromosome.genes[param_numbers[i]] * rad2deg;
				file << std::endl;
			}	
			else {
				file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << chromosome.genes[param_numbers[i]];
				file << std::endl;
			}
			if ((i == 27) && (param_modes[27] == 1)) {
				for (std::vector<int>::size_type j = param_numbers[i] + 1; j != chromosome.genes.size(); j++) {
					file << std::setw(col_width) << " ";
					file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << chromosome.genes[j];
					file << std::endl;
				}
			}
		}
	}
	file.close();
}

void GeneticAlgorithm::record_symmetric_parameters(PeldorCalculator const& peldor_calc,
	                                               Chromosome const& chromosome,
	                                               std::vector<exp_signal> const& signals_exp,
												   Spin const& spinA, Spin const& spinB,
												   std::vector<int> const& param_numbers,
												   std::vector<int> const& param_modes,
												   std::vector<bound> const& param_bounds,
												   genetic_parameters const& genetic_param,
												   output_parameters const& output_param) const
{

	// Read out the optimized values of angles
	double xi(0), phi(0), alpha(0), betta(0), gamma(0);
	if (param_numbers[2] != -1) xi = chromosome.genes[param_numbers[2]];
	if (param_numbers[4] != -1) phi = chromosome.genes[param_numbers[4]];
	if (param_numbers[6] != -1) alpha = chromosome.genes[param_numbers[6]];
	if (param_numbers[8] != -1) betta = chromosome.genes[param_numbers[8]];
	if (param_numbers[10] != -1) gamma = chromosome.genes[param_numbers[10]];
	// Create the array for the symmetry-related sets of parameters
	const size_t n_sym = 16;
	std::vector<double> xis; xis.reserve(n_sym);
	std::vector<double> phis; phis.reserve(n_sym);
	std::vector<double> alphas; alphas.reserve(n_sym);
	std::vector<double> bettas; bettas.reserve(n_sym);
	std::vector<double> gammas; gammas.reserve(n_sym);
	xis.push_back(xi);
	phis.push_back(phi);
	alphas.push_back(alpha);
	bettas.push_back(betta);
	gammas.push_back(gamma);
	// Calculate the orientation of the distance vector
	std::vector<double> dist_vec; dist_vec.reserve(3);
	dist_vec.push_back( cos(phi)*sin(xi) );
	dist_vec.push_back( sin(phi)*sin(xi) );
	dist_vec.push_back( cos(xi) );
	// Calculate the rotation matrix from frame A to frame B
	ZXZEulerRM const* R = new ZXZEulerRM(alpha,betta,gamma);
	// Create the 180° rotation matrices
	ZXZEulerRM const* Rx = new ZXZEulerRM(0.0, PI, 0.0);
	ZXZEulerRM const* Ry = new ZXZEulerRM(PI, PI, 0.0);
	ZXZEulerRM const* Rz = new ZXZEulerRM(PI, 0.0, 0.0);
	// Calculate the symmetry-related sets of parameters
	std::vector<double> dist_vec_sym; dist_vec_sym.reserve(3);
	double xi_sym(0), phi_sym(0), alpha_sym(0), betta_sym(0), gamma_sym(0);
	double R_sym[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
	// Rotation of spin A about gx
	//----------------------------
	dist_vec_sym = Rx->rotate_vector(dist_vec);
	xi_sym = acos(dist_vec_sym[2]);
    if ((xi_sym == 0) || (xi_sym == PI)) {
		phi_sym = 0;
	}
    else {
        if ( dist_vec_sym[0] == 0 ) {
            if ( dist_vec_sym[1] == sin(xi_sym) ) phi_sym = 0.5*PI;
            else                                  phi_sym = -0.5*PI;
		}
        else phi_sym = atan2(dist_vec_sym[1], dist_vec_sym[0]);
	}
	if (phi_sym < 0) phi_sym += 2*PI;
	xis.push_back(xi_sym);
	phis.push_back(phi_sym);
	alphas.push_back(alpha);
	bettas.push_back(betta);
	gammas.push_back(gamma);
	dist_vec_sym.clear();
	// Rotation of spin A about gy
	//----------------------------
	dist_vec_sym = Ry->rotate_vector(dist_vec);
	xi_sym = acos(dist_vec_sym[2]);
    if ((xi_sym == 0) || (xi_sym == PI)) {
		phi_sym = 0;
	}
    else {
        if ( dist_vec_sym[0] == 0 ) {
            if ( dist_vec_sym[1] == sin(xi_sym) ) phi_sym = 0.5*PI;
            else                                  phi_sym = -0.5*PI;
		}
        else phi_sym = atan2( dist_vec_sym[1], dist_vec_sym[0] );
	}
	if (phi_sym < 0) phi_sym += 2*PI;
	xis.push_back(xi_sym);
	phis.push_back(phi_sym);
	alphas.push_back(alpha);
	bettas.push_back(betta);
	gammas.push_back(gamma);
	dist_vec_sym.clear();
	// Rotation of spin A about gz
	//----------------------------
	dist_vec_sym = Rz->rotate_vector(dist_vec);
	xi_sym = acos(dist_vec_sym[2]);
    if ((xi_sym == 0) || (xi_sym == PI)) {
		phi_sym = 0;
	}
    else {
        if ( dist_vec_sym[0] == 0 ) {
            if ( dist_vec_sym[1] == sin(xi_sym) ) phi_sym = 0.5*PI;
            else                                  phi_sym = -0.5*PI;
		}
        else phi_sym = atan2( dist_vec_sym[1], dist_vec_sym[0] );
	}
	if (phi_sym < 0) phi_sym += + 2*PI;
	xis.push_back(xi_sym);
	phis.push_back(phi_sym);
	alphas.push_back(alpha);
	bettas.push_back(betta);
	gammas.push_back(gamma);
	dist_vec_sym.clear();
	// Rotation of spin B about gx
	//----------------------------
	R->multiply_by_matrix(Rx->R, R_sym);
	betta_sym = acos(R_sym[2][2]);
    if (betta_sym == 0) {
        gamma_sym = 0;
        alpha_sym = atan2(R_sym[1][0], R_sym[1][1]);
	}
    else {
        gamma_sym = atan2(R_sym[2][0]/sin(betta_sym), R_sym[2][1]/sin(betta_sym));
        alpha_sym = atan2(R_sym[0][2]/sin(betta_sym), -R_sym[1][2]/sin(betta_sym));
	}
    if (alpha_sym < 0) alpha_sym += 2*PI;
    if (gamma_sym < 0) gamma_sym += 2*PI;
	xis.push_back(xi);
	phis.push_back(phi);
	alphas.push_back(alpha_sym);
	bettas.push_back(betta_sym);
	gammas.push_back(gamma_sym);
	// Rotation of spin B about gy
	//----------------------------
	R->multiply_by_matrix(Ry->R, R_sym);
	betta_sym = acos(R_sym[2][2]);
    if (betta_sym == 0) {
        gamma_sym = 0;
        alpha_sym = atan2(R_sym[1][0], R_sym[1][1]);
	}
    else {
        gamma_sym = atan2(R_sym[2][0]/sin(betta_sym), R_sym[2][1]/sin(betta_sym));
        alpha_sym = atan2(R_sym[0][2]/sin(betta_sym), -R_sym[1][2]/sin(betta_sym));
	}
    if (alpha_sym < 0) alpha_sym += 2*PI;
    if (gamma_sym < 0) gamma_sym += 2*PI;
	xis.push_back(xi);
	phis.push_back(phi);
	alphas.push_back(alpha_sym);
	bettas.push_back(betta_sym);
	gammas.push_back(gamma_sym);
	// Rotation of spin B about gz
	//----------------------------
	R->multiply_by_matrix(Rz->R, R_sym);
	betta_sym = acos(R_sym[2][2]);
    if (betta_sym == 0) {
        gamma_sym = 0;
        alpha_sym = atan2(R_sym[1][0], R_sym[1][1]);
	}
    else {
        gamma_sym = atan2(R_sym[2][0]/sin(betta_sym), R_sym[2][1]/sin(betta_sym));
        alpha_sym = atan2(R_sym[0][2]/sin(betta_sym), -R_sym[1][2]/sin(betta_sym));
	}
    if (alpha_sym < 0) alpha_sym += 2*PI;
    if (gamma_sym < 0) gamma_sym += 2*PI;
	xis.push_back(xi);
	phis.push_back(phi);
	alphas.push_back(alpha_sym);
	bettas.push_back(betta_sym);
	gammas.push_back(gamma_sym);
	// Combined rotations
	//----------------------------
	for (size_t i = 1; i < 4; ++i) {
		for (size_t j = 4; j < 7; ++j) {
			xis.push_back(xis[i]);
			phis.push_back(phis[i]);
			alphas.push_back(alphas[j]);
			bettas.push_back(bettas[j]);
			gammas.push_back(gammas[j]);
		}
	}
	// Store the calculated symmetry-related sets of parameters in a vector
	std::vector<std::vector<double>> genes_sym; genes_sym.reserve(n_sym);
	for (size_t i = 0; i < n_sym; ++i) {
		genes_sym.push_back(chromosome.genes);
		if (param_numbers[2] != -1)  genes_sym[i][param_numbers[2]] = xis[i];
		if (param_numbers[4] != -1)  genes_sym[i][param_numbers[4]] = phis[i];
		if (param_numbers[6] != -1)  genes_sym[i][param_numbers[6]] = alphas[i];
		if (param_numbers[8] != -1)  genes_sym[i][param_numbers[8]] = bettas[i];
		if (param_numbers[10] != -1) genes_sym[i][param_numbers[10]] = gammas[i];
	}
	// Calculate the score for the symmetry-related sets of parameters
	genetic_parameters genetic_param_sym(genetic_param);
	genetic_param_sym.size_generation = n_sym;
	Generation* generation_sym = new Generation(genetic_param_sym);
	for (size_t i = 0; i < n_sym; ++i) {
		generation_sym->chromosomes.push_back(chromosome);
		generation_sym->chromosomes[i].genes = genes_sym[i];
	}
	generation_sym->score_chromosomes(peldor_calc, signals_exp, spinA, spinB, param_numbers, param_modes, genetic_param);
	// Create & open a file
	std::ostringstream filename;
	filename << output_param.directory << "symmetric_parameters.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the names of columns
	size_t const col_width1 = 30;
	size_t const col_width2 = 10;
	file << std::left << std::setw(col_width1) << "Parameter";
	for (size_t i = 0; i < n_sym; ++i) file << std::left << std::setw(col_width2) << "Value";
	file << std::endl;
	// Fill in the columns with names of parameters and their values
	const char *param_names[] = { "distance mean 1 (nm)","distance width 1 (nm)","xi mean 1 (deg)","xi width 1 (deg)","phi mean 1 (deg)","phi width 1 (deg)",
		                          "alpha mean 1 (deg)","alpha width 1 (deg)","betta mean 1 (deg)","betta width 1 (deg)","gamma mean 1 (deg)","gamma width 1 (deg)",
								  "distance mean 2 (nm)", "distance width 2 (nm)", "xi mean 2 (deg)", "xi width 2 (deg)", "phi mean 2 (deg)", "phi width 2 (deg)",
								  "alpha mean 2 (deg)", "alpha width 2 (deg)", "betta mean 2 (deg)", "betta width 2 (deg)", "gamma mean 2 (deg)", "gamma width 2 (deg)",
								  "ratio 1 to 2", "J mean (MHz)", "J width (MHz)", "inv. efficiency" };
	for (size_t i = 0; i < 28; ++i) {
		file << std::left << std::setw(col_width1) << param_names[i];
		if (param_numbers[i] == -1) {
			file << std::setw(col_width2*n_sym) << " ";
			file << std::endl;
		}
		else {
			if (std::count(angular_indices, angular_indices + 20, i)) {
				for (size_t k = 0; k < n_sym; ++k) {
					file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(0) << generation_sym->chromosomes[k].genes[param_numbers[i]] * rad2deg;
				}
				file << std::endl;
			}	
			else {
				for (size_t k = 0; k < n_sym; ++k) {
					file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(2) << generation_sym->chromosomes[k].genes[param_numbers[i]];
				}
				file << std::endl;
			}
			if ((i == 27) && (param_modes[27] == 1)) {
				for (std::vector<int>::size_type j = param_numbers[i] + 1; j != chromosome.genes.size(); j++) {
					file << std::setw(col_width1) << " ";
					for (size_t k = 0; k < n_sym; ++k) {
						file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(2) << generation_sym->chromosomes[k].genes[j];
					}
					file << std::endl;
				}
			}
		}
		
	}
	if (genetic_param.merit_function == 1)      file << std::left << std::setw(col_width1) << "RMSD";
	else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width1) << "RMSD/PCC";
	else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width1) << "PCC";
	for (size_t k = 0; k < n_sym; ++k) {
		file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(4) << generation_sym->chromosomes[k].fitness;
	}
	file.close();
	delete R;
	delete Rx;
	delete Ry;
	delete Rz;
	delete generation_sym;
}

void GeneticAlgorithm::record_form_factor(PeldorCalculator const& peldor_calc,
	                                      Chromosome const& chromosome,
										  std::vector<exp_signal> const& signals_exp,
										  Spin const& spinA, Spin const& spinB,
										  std::vector<int> const& param_numbers,
										  std::vector<int> const& param_modes,
										  output_parameters const& output_param) const
{
	const size_t numOfSignals = signals_exp.size();
	std::vector<double> form_factor; form_factor.reserve(90);
	std::vector<std::vector<double>> form_factors; form_factors.reserve(numOfSignals);
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 engine(seq);
	for (size_t j = 0; j < numOfSignals; ++j) {
		// Calculate the form-factor for the best chromosome and different frequency offsets
		form_factor = peldor_calc.calculate_form_factor(chromosome.genes, signals_exp[j], spinA, spinB, 
				                                        param_numbers, param_modes, engine, j);
		form_factors.push_back(form_factor);
		form_factor.clear();
	}
	// Create & open a file
	std::ostringstream filename;
	filename << output_param.directory << "form_factor.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the names of columns
	size_t const col_width = 15;
	file << std::left << std::setw(col_width) << "Angle(deg)";
	std::ostringstream ff_name;
	for (size_t i = 0; i < signals_exp.size(); ++i) {
		ff_name << "Form-factor" << i + 1;
		file << std::left << std::setw(col_width) << ff_name.str();
		ff_name.str(std::string()); ff_name.clear();
	}
	file << std::endl;
	// Write the calculated form factors
	for (size_t t = 0; t < 90; ++t) {
		file << std::left << std::setw(col_width) << std::fixed << std::setprecision(1) << static_cast<double>(t) + 0.5;
		for (size_t i = 0; i < form_factors.size(); ++i) {
			file << std::left << std::setw(col_width) << std::fixed << std::setprecision(5) << form_factors[i][t];
		}
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_error_plot(PeldorCalculator const& peldor_calc,
	                                     Chromosome const& chromosome,
									     std::vector<exp_signal> const& signals_exp,
									     Spin const& spinA, Spin const& spinB,
									     std::vector<int> const& param_numbers,
									     std::vector<int> const& param_modes,
										 std::vector<bound> const& param_bounds,
										 genetic_parameters const& genetic_param,
									     output_parameters const& output_param) const
{
	const size_t n_plots = output_param.error_plot_variables.size();
	const size_t plot_size = output_param.error_plot_size;
	size_t n_vars(0);
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 engine(seq);
	for (size_t i = 0; i < n_plots; ++i) { // iterate over the error plots
		n_vars = output_param.error_plot_variables[i].size();
		// Check if the selected variables of the error plot were optimized
		int indexVariable(0);
		bool inappropriateVariables(0);
		for (size_t j = 0; j < n_vars; ++j) {
			indexVariable = output_param.error_plot_variables[i][j] - 1;
			if (param_numbers[indexVariable] == -1) {
				std::cout << "  Inappropriate variable (it was not optimized): " << output_param.error_plot_variables[i][j] << std::endl;
				inappropriateVariables = 1;
				break;
			}
		}
		if (inappropriateVariables) continue;
		// Create a generation with different values of the parameters to be varied
		genetic_parameters genetic_param_ep(genetic_param);
		genetic_param_ep.size_generation = output_param.error_plot_size;
		Generation* generation_ep = new Generation(genetic_param_ep);
		for (size_t k = 0; k < plot_size; ++k) generation_ep->chromosomes.push_back(chromosome); // initialize all chromosomes with the optimized chromosome
		int indexGene(0);
		for (size_t j = 0; j < n_vars; ++j) { // vary the genes corresponding to the defined variables
			indexGene = param_numbers[output_param.error_plot_variables[i][j] - 1];
			for (size_t k = 0; k < plot_size; ++k) {
				generation_ep->chromosomes[k].genes[indexGene] = chromosome.create_random_gene(param_bounds[indexGene].lower, param_bounds[indexGene].upper, engine);
			}
		}
		// Score the generation
		generation_ep->score_chromosomes(peldor_calc, signals_exp, spinA, spinB, param_numbers, param_modes, genetic_param);
		generation_ep->sort_chromosomes();
		// Create & open a file
		std::ostringstream filename;
		filename << output_param.directory << "error_plot";
		for (size_t j = 0; j < n_vars; ++j) filename << "_" << output_param.error_plot_variables[i][j];
		filename << ".dat";
		std::fstream file;
		file.open(filename.str().c_str(), std::ios::out);
		// Write the names of columns
		size_t const col_width = 15;
		std::ostringstream col_name;
		for (size_t j = 0; j < n_vars; ++j) {
			col_name << "Parameter" << (j + 1);
			file << std::left << std::setw(col_width) << col_name.str();
			col_name.str(std::string()); col_name.clear();
		}
		if (genetic_param.merit_function == 1)      file << std::left << std::setw(col_width) << "RMSD" << std::endl;
		else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width) << "RMSD/PCC" << std::endl;
		else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width) << "PCC" << std::endl;
		// Write the error plot
		for (size_t k = 0; k < plot_size; ++k) {
			for (size_t j = 0; j < n_vars; ++j) { 
				indexGene = param_numbers[output_param.error_plot_variables[i][j] - 1];
				if (std::count(angular_indices, angular_indices + 20, indexGene)) {
					file << std::left << std::setw(col_width) << std::fixed << std::setprecision(0) << generation_ep->chromosomes[k].genes[indexGene] * rad2deg;
				}
				else {
					file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << generation_ep->chromosomes[k].genes[indexGene];
				}
			}
			file << std::left << std::setw(col_width) << std::fixed << std::setprecision(5) << generation_ep->chromosomes[k].fitness;
			file << std::endl;
		}
		file.close();
		// Clean up
		delete generation_ep;
	}
}

void GeneticAlgorithm::run_error_plot(std::vector<exp_signal> const& signals_exp,
	                                  Spin const& spinA, Spin const& spinB,
									  std::vector<int> const& param_numbers,
									  std::vector<int> const& param_modes,
									  std::vector<bound> const& param_bounds,
									  genetic_parameters& genetic_param,
									  output_parameters& output_param,
									  errorplot_parameters const& errorplot_param) const
{
	// Create a chromosome and set its genes to the read values
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 engine(seq);
	Chromosome* chromosome = new Chromosome(param_bounds,engine);
	chromosome->genes = errorplot_param.optimized_genes;
	// Initialize the calculator of PELDOR signals
	PeldorCalculator* peldor_calc = new PeldorCalculator(genetic_param.num_avg);
	peldor_calc->calculate_spinA_excitation(signals_exp, spinA);
	// Re-define the output parameters
	output_param.directory = errorplot_param.output_directory;
	output_param.error_plot_variables = errorplot_param.error_plot_variables;
	output_param.error_plot_size = errorplot_param.error_plot_size;
	// Record the error plot
	record_error_plot(*peldor_calc, *chromosome, signals_exp, spinA, spinB, param_numbers, param_modes, 
		              param_bounds, genetic_param, output_param);
	// Clean up
	delete peldor_calc;
	delete chromosome;
}