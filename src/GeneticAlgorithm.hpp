#pragma once

#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include "Utilities.hpp"
#include "Individual.hpp"
#include <memory>
#include <vector>
#include <random>

using CrossoverFunction = void(*)(
    std::vector<std::unique_ptr<Individual>> &,
    std::mt19937 &,
    const int,
    const int,
    const int
);

using MutationFunction = void(*)(
    const std::vector <std::unique_ptr <Individual>> &,
    const std::vector <size_t> &,
    const double,
    const double,
    const double,
    const double,
    std::mt19937 &
);

std::unique_ptr <std::vector <double>> generate_solution(
    const Func & function,
    const int multiplier,
    const int constant,
    const double step
);

std::unique_ptr <Individual> generate_solution(
    std::mt19937 & generator,
    const double min_value,
    const double max_value,
    const double step
);

void crossover_type1(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
);

void crossover_type2(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
);

void crossover_type3(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
);

void crossover_type4(
    const Individual & parent1,
    const Individual & parent2,
    std::vector<double> & child1_solution,
    std::vector<double> & child2_solution,
    std::mt19937 & generator
);

std::pair <std::unique_ptr <Individual>, std::unique_ptr <Individual>> execute_crossover(
    const Individual & parent1,
    const Individual & parent2,
    std::mt19937 & generator,
    const int crossover_strategy
);

static void execute_parallel_crossover(
    std::vector <std::unique_ptr <Individual>> & population,
    const int current_population_size,
    const int crossover_strategy,
    const size_t index
);

void parallel_crossover(
    std::vector <std::unique_ptr <Individual>> & population,
    std::mt19937 & generator,
    const int population_size,
    const int max_population_size,
    const int crossover_strategy
);

void single_thread_crossover(
    std::vector <std::unique_ptr<Individual>> & population,
    std::mt19937 & generator,
    const int population_size,
    const int max_population_size,
    const int crossover_strategy
);

std::pair <int, int> choose_crossover_candidates(
    std::mt19937 & generator,
    const int population_size
);

std::unique_ptr <std::vector <size_t>> choose_mutation_candidates(
    std::mt19937 & generator,
    const size_t max_population_size,
    const double mutation_rate
);

void execute_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const double average_fitness_value,
    const double min_range,
    const double max_range,
    const double step,
    std::mt19937 & generator
);

void parallel_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const double average_fitness_value,
    const double min_range,
    const double max_range,
    const double step,
    std::mt19937 & generator
);

double evaluate_average_fitness(
    const std::vector <std::unique_ptr <Individual>> & population
);

double evaluate_maximum_fitness_value(
    const std::vector <std::unique_ptr <Individual>> & population
);

void print_cycle_data(
    const int cycle,
    const double average_fitness_value,
    const double max_fitness_value,
    const int mutation_amount
);

void ga_loop(
    const std::unique_ptr <std::vector <std::unique_ptr <Individual>>> & population,
    const std::unique_ptr <std::vector <double>> & target,
    const Config & config,
    std::mt19937 & generator
);

#endif
