#include "Config.hpp"
#include "Utilities.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

Config parse_config(
    const std::string & file_path
) {
    Config config;
    boost::property_tree::ptree ptree;
    boost::property_tree::ini_parser::read_ini(file_path, ptree);

    config.step = ptree.get <double> ("step");
    config.population_size = ptree.get <int> ("population_size");
    config.max_population_size = ptree.get <int> ("max_population_size");
    config.mutation_rate = ptree.get <double> ("mutation_rate");
    config.cycles = ptree.get <int> ("cycles");
    config.crossover_strategy = ptree.get <int> ("crossover_strategy");
    config.constant = ptree.get <int> ("function.constant");
    config.function = ptree.get <std::string> ("function.function");
    config.argument_multiplier = ptree.get <int> ("function.argument_multiplier");
    config.min_range = ptree.get <double> ("min_range");
    config.max_range = ptree.get <double> ("max_range");

    validate_config_step(config.step);
    validate_population_size(config.population_size, config.max_population_size);
    validate_mutation_rate(config.mutation_rate);
    validate_cycles(config.cycles);
    validate_crossover_strategy(config.crossover_strategy);
    validate_function(config.constant, config.function, config.argument_multiplier);
    validate_range(config.min_range, config.max_range);

    return config;
}

void print_config(
    const Config & config
) {
    std::cout << "step: " << config.step << '\n';
    std::cout << "population size: " << config.population_size << '\n';
    std::cout << "max population size: " << config.max_population_size << '\n';
    std::cout << "mutation rate: " << config.mutation_rate << '\n';
    std::cout << "cycles: " << config.cycles << '\n';
    std::cout << "crossover strategy: " << config.crossover_strategy << '\n';
    std::cout << "function: ";
    if (config.constant != 0) {
        std::cout << config.constant << " + ";
    }

    if (config.function == "square") {
        std::cout << config.function << "(x) * " << config.argument_multiplier << '\n';
    }
    else {
        std::cout << config.function << "(";
        if (config.argument_multiplier != 1) {
            std::cout << config.argument_multiplier << "x)\n";
        }
        else {
            std::cout << "x)\n";
        }
    }
}