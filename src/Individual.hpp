#pragma once

#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include "Config.hpp"
#include <memory>
#include <string>
#include <vector>
#include <random>

struct Individual {
    double fitness;
    std::unique_ptr <std::vector <double>> solution;

    Individual(double fit, std::unique_ptr <std::vector <double>> sol) : fitness(fit), solution(std::move(sol)) {}
};

std::unique_ptr <std::vector <std::unique_ptr <Individual>>> create_population(
    const Config & config,
    std::mt19937 & generator
);

void evaluate_fitness(
    const std::vector <double> & target,
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