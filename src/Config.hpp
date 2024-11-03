#pragma once

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

struct Config {
    double step;
    int population_size;
    int max_population_size;
    double mutation_rate;
    int cycles;
    int crossover_strategy;
    int constant;
    std::string function;
    int argument_multiplier;
    double min_range;
    double max_range;
    int computation_mode;
};

Config parse_config(
    const std::string & file_path
);

void print_config(
    const Config & config
);

#endif
