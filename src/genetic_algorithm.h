#ifndef GENETIC_ALGORITH_H
#define GENETIC_ALGORITH_H

#include "definitions.h"
#include "spin.h"
#include "peldor_calculator.h"
#include <vector>
#include <cmath>
#include <random>

class Chromosome
{
public:
	size_t size;
	std::vector<double> genes;
	double fitness;

	Chromosome(std::vector<bound> const& param_bounds,
			   std::mt19937& urng);

	// Create a random gene
	double create_random_gene(double const& gene_min,
		                      double const& gene_max,
							  std::mt19937& urng) const;

	// Overload less operator
	bool operator < (const Chromosome& A) const { return this->fitness < A.fitness; };
};

class Generation
{
public:
	size_t size;
	std::vector<Chromosome> chromosomes;
	std::vector<Chromosome> offspring;

	Generation(genetic_parameters const& genetic_param);

	// Create the first generation
	void create_initial_generation(std::vector<bound> const& param_bounds,
		                           genetic_parameters const& genetic_param);

	// Calculate the fitness of individual genes
	void score_chromosomes(PeldorCalculator const& peldor_calc,
		                   std::vector<exp_signal> const& signals_exp,
						   Spin const& spinA, Spin const& spinB,
						   std::vector<int> const& param_numbers,
						   std::vector<int> const& param_modes,
						   genetic_parameters const& genetic_param);

	// Sort chromosomes based on the fitness of the genes
	void sort_chromosomes();

	// Create new chromosomes from the old chromosomes
	void produce_offspring(genetic_parameters const& genetic_param, std::vector<bound> const& param_bounds);

	int tournament_selection(std::mt19937& urng) const;

	void crossover_chromosomes(Chromosome& chromosome1, Chromosome& chromosome2, 
		                       double const& prob_crossover,
							   std::mt19937& urng) const;

	void mutate_chromosome(Chromosome& chromosome, 
		                   double const& prob_mutation,
		                   std::vector<bound> const& param_bounds, 
						   std::mt19937& urng) const;
};

class GeneticAlgorithm
{
public:
	GeneticAlgorithm();

	void run_optimization(std::vector<exp_signal> const& signals_exp,
	                      Spin const& spinA, Spin const& spinB,
						  std::vector<int> const& param_numbers,
						  std::vector<int> const& param_modes,
						  std::vector<bound> const& param_bounds,
						  genetic_parameters const& genetic_param,
						  output_parameters const& output_param) const;

	// Save the fitness of the best chromosome
	void record_score(std::vector<double> const& best_fitness,
		              genetic_parameters const& genetic_param,
					  output_parameters const& output_param) const;

	// Save the calculated fit to the PELDOR signals for the best chromosome
	void record_fit(PeldorCalculator const& peldor_calc,
		            Chromosome const& chromosome,
		            std::vector<exp_signal> const& signals_exp,
					Spin const& spinA, Spin const& spinB,
					std::vector<int> const& param_numbers,
					std::vector<int> const& param_modes,
					output_parameters const& output_param) const;

	// Save the genes of the best chromosome
	void record_parameters(Chromosome const& chromosome,
		                   std::vector<int> const& param_numbers,
						   std::vector<int> const& param_modes,
						   output_parameters const& output_param) const;

	// Save the symmetry-related sets of fitting parameters
	void record_symmetric_parameters(PeldorCalculator const& peldor_calc,
		                             Chromosome const& chromosome,
									 std::vector<exp_signal> const& signals_exp,
									 Spin const& spinA, Spin const& spinB,
									 std::vector<int> const& param_numbers,
									 std::vector<int> const& param_modes,
									 std::vector<bound> const& param_bounds,
									 genetic_parameters const& genetic_param,
									 output_parameters const& output_param) const;

	// Save the calculated form factors for the PELDOR signals
	void record_form_factor(PeldorCalculator const& peldor_calc,
		                    Chromosome const& chromosome,
							std::vector<exp_signal> const& signals_exp,
							Spin const& spinA, Spin const& spinB,
							std::vector<int> const& param_numbers,
							std::vector<int> const& param_modes,
							output_parameters const& output_param) const;
	
	// Record the error plot
	void record_error_plot(PeldorCalculator const& peldor_calc,
		                   Chromosome const& chromosome,
						   std::vector<exp_signal> const& signals_exp,
						   Spin const& spinA, Spin const& spinB,
						   std::vector<int> const& param_numbers,
						   std::vector<int> const& param_modes,
						   std::vector<bound> const& param_bounds,
						   genetic_parameters const& genetic_param,
						   output_parameters const& output_param) const;

	// Record only error plot using previous fitting results
	void run_error_plot(std::vector<exp_signal> const& signals_exp,
		                Spin const& spinA, Spin const& spinB,
						std::vector<int> const& param_numbers,
						std::vector<int> const& param_modes,
						std::vector<bound> const& param_bounds,
						genetic_parameters& genetic_param,
						output_parameters& output_param,
						errorplot_parameters const& errorplot_param) const;
};

#endif