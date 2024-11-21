#include "Utilities.hpp"
#include <cmath>
#include <set>
#include <iostream>
#include <stdexcept>

std::map <std::string, Func> string_to_function_map = {
    {"sin", std::sin},
    {"cos", std::cos},
    {"square", square},
    {"ln", log_zero_allowed},
    {"sqrt", std::sqrt},
};

float square(const float x) {
    return std::pow(x, 2);
}

float log_zero_allowed(const float x) {
    if (x == 0.0) {
        return std::log(0.001);
    }

    return std::log(x);
}

void validate_config_step(float step) {
    std::set <float> valid_steps = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5};

    if (valid_steps.find(step) == valid_steps.end()) {
        std::cerr << "invalid step value\n";
        throw std::invalid_argument("invalid step value");
    }
}

void validate_population_size(int population_size, int max_population_size) {
    int minimal_max_size = population_size + 2;

    if (max_population_size < minimal_max_size) {
        std::cerr << "invalid population and / or max population size value\n";
        throw std::invalid_argument("invalid population and / or max population size value\n");
    }

    if ((population_size % 2 == 0) && (max_population_size % 2 != 0)) {
        std::cerr << "invalid population and / or max population size value\n";
        throw std::invalid_argument("invalid population and / or max population size value\n");
    }

    if ((population_size % 2 != 0) && (max_population_size % 2 == 0)) {
        std::cerr << "invalid population and / or max population size value\n";
        throw std::invalid_argument("invalid population and / or max population size value\n");
    }
}

void validate_mutation_rate(float mutation_rate) {
    if (mutation_rate < 0.0 || mutation_rate > 1.0) {
        std::cerr << "invalid mutation rate value\n";
        throw std::invalid_argument("invalid mutation rate value\n");
    }
}

void validate_cycles(int cycles) {
    if (cycles < 0) {
        std::cerr << "invalid cycles value\n";
        throw std::invalid_argument("invalid cycles value\n");
    }
}

void validate_crossover_strategy(int crossover_strategy) {
    if (crossover_strategy < 0 || crossover_strategy > 3) {
        std::cerr << "invalid crossover strategy value\n";
        throw std::invalid_argument("invalid crossover strategy value\n");
    }
}

void validate_function(int constant, std::string function, int argument_multiplier) {
    if (constant != 0 && constant != 1) {
        std::cerr << "invalid function constant value\n";
        throw std::invalid_argument("invalid function constant value\n");
    }

    auto string_to_function_map_iterator = string_to_function_map.find(function);
    if (string_to_function_map_iterator == string_to_function_map.end()) {
        std::cerr << "invalid function name\n";
        throw std::invalid_argument("invalid function name\n");
    }

    if (argument_multiplier < 1 || argument_multiplier > 4) {
        std::cerr << "invalid function argument multiplier value\n";
        throw std::invalid_argument("invalid function argument multiplier value\n");
    }
}

void validate_range(float min_range, float max_range) {
    if (min_range >= max_range) {
        std::cerr << "invalid range values\n";
        throw std::invalid_argument("invalid range values\n");
    }
}

void validate_computation_mode(
    int computation_mode
) {
    if ((computation_mode < 0) || (computation_mode > 2)) {
        std::cerr << "invalid computation mode\n";
        throw std::invalid_argument("invalid computation mode\n");
    }
}

void validate_print_interval(
    int print_interval
) {
    if (print_interval <= 0) {
        std::cerr << "invalid print interval\n";
        throw std::invalid_argument("invalid print interval\n");
    }
}

void validate_approximation_tolerance(
    float approximation_tolerance
) {
    if (approximation_tolerance <= 0.0) {
        std::cerr << "invalid approximation tolerance\n";
        throw std::invalid_argument("invalid approximation tolerance\n");
    }
}
