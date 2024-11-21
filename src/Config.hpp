#pragma once

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

struct Config {
    double step;
    int population_size;
    int max_population_size;
    float mutation_rate;
    int cycles;
    int crossover_strategy;
    int constant;
    std::string function;
    int argument_multiplier;
    float min_range;
    float max_range;
    int computation_mode;
    int print_interval;
    float approximation_tolerance;
};

Config parse_config(
    const std::string & file_path
);

void print_config(
    const Config & config
);

#endif
