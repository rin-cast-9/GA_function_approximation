#pragma once

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>
#include <map>

using Func = float (*)(float);

extern std::map <std::string, Func> string_to_function_map;

float square(
    const float x
);

float log_zero_allowed(
    const float x
);

void validate_config_step(
    float step
);

void validate_population_size(
    int population_size,
    int max_population_size
);

void validate_mutation_rate(
    float mutation_rate
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
    float min_range,
    float max_range
);

void validate_computation_mode(
    int computation_mode
);

void validate_print_interval(
    int print_interval
);

void validate_approximation_tolerance(
    float approximation_tolerance
);

#endif
