#pragma once

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>
#include <map>
#include <cmath>

using Func = double (*)(double);

extern std::map <std::string, Func> string_to_function_map;

double square(
    const double x
);

void validate_config_step(
    double step
);

void validate_population_size(
    int population_size,
    int max_population_size
);

void validate_mutation_rate(
    double mutation_rate
);

void validate_cycles(
    int cycles
);

void validate_crossover_strategy(
    int crossover_strategy
);

void validate_function(
    int constant,
    std::string function,
    int argument_multiplier
);

void validate_range(
    double min_range,
    double max_range
);

#endif