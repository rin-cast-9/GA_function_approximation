#pragma once

#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include "Config.hpp"
#include <memory>
#include <vector>
#include <random>

struct Individual {
    float fitness;
    std::unique_ptr <std::vector <float>> solution;

    Individual(float fit, std::unique_ptr <std::vector <float>> sol) : fitness(fit), solution(std::move(sol)) {}
};

using FitnessFunction = void(*)(
    const std::vector <float> &,
    const std::vector <std::unique_ptr <Individual>> &
);

std::unique_ptr <std::vector <std::unique_ptr <Individual>>> create_population(
    const Config & config,
    std::mt19937 & generator
);

void evaluate_fitness(
    const std::vector <float> & target,
    const std::vector <std::unique_ptr <Individual>> & population
);

void parallel_evaluate_fitness(
    const std::vector <float> & target,
    const std::vector <std::unique_ptr <Individual>> & population
);

void print_population(
    const std::vector <std::unique_ptr <Individual>> & population
);

void perform_selection(
    std::vector <std::unique_ptr <Individual>> & population,
    const int population_size
);

#endif
