#include "Individual.hpp"
#include "GeneticAlgorithm.hpp"
#include <iostream>
#include <algorithm>

std::unique_ptr <std::vector <std::unique_ptr <Individual>>> create_population(
    const Config & config, std::mt19937 & generator
) {
    auto population = std::make_unique <std::vector <std::unique_ptr <Individual>>> ();

    for (int i = 0; i < config.population_size; ++ i) {
        population->push_back(generate_solution(generator, config.min_range, config.max_range, config.step));
    }

    return population;
}

void evaluate_fitness(
    const std::vector <double> & target,
    const std::vector <std::unique_ptr <Individual>> & population
) {
    size_t population_size = population.size();
    size_t individual_size = target.size();

    auto fitness_values = std::make_unique <std::vector <double>> (population_size);

    for (int i = 0; i < population_size; ++ i) {
        double individual_fitness_value = 0.0;
        auto & individuals_solution = * population[i]->solution;
        for (int j = 0; j < individual_size; ++ j) {
            individual_fitness_value += std::abs(target[j] - individuals_solution[j]);
        }
        population[i]->fitness = individual_fitness_value;
    }
}

void print_population(
    const std::vector <std::unique_ptr <Individual>> & population
) {
    for (const auto & individual : population) {
        std::cout << "[" << individual->fitness << "] ";

        for (const auto & element : * individual->solution) {
            std::cout << element << ' ';
        }
    }
    std::cout << "]\n";
}

void perform_selection(
    std::vector <std::unique_ptr <Individual>> & population,
    const int population_size    
) {
    std::sort(
        population.begin(),
        population.end(),
        [](const std::unique_ptr <Individual> & a, const std::unique_ptr <Individual> & b) {
            return a->fitness < b->fitness;
        }
    );

    population.resize(population_size);
}