#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


double square(const double x) {
    return std::pow(x, 2);
}

using Func = double (*)(double);

std::map <std::string, Func> string_to_function_map = {
    {"sin", std::sin},
    {"cos", std::cos},
    {"square", square},
    {"ln", std::log},
    {"sqrt", std::sqrt},
};

struct Config {
    double step;
    int population_size;
    double mutation_rate;
    int cycles;
    int constant;
    std::string function;
    int argument_multiplier;
};

struct Individual {
    double fitness;
    std::unique_ptr <std::vector <double>> solution;

    Individual(double fit, std::unique_ptr <std::vector <double>> sol) : fitness(fit), solution(std::move(sol)) {}
};

Config parse_config(const std::string & file_path) {
    Config config;
    boost::property_tree::ptree ptree;
    boost::property_tree::ini_parser::read_ini(file_path, ptree);

    config.step = ptree.get <double> ("step");
    config.population_size = ptree.get <int> ("population_size");
    config.mutation_rate = ptree.get <double> ("mutation_rate");
    config.cycles = ptree.get <int> ("cycles");
    config.constant = ptree.get <int> ("function.constant");
    config.function = ptree.get <std::string> ("function.function");
    config.argument_multiplier = ptree.get <int> ("function.argument_multiplier");

    return config;
}

void print_config(const Config & config) {
    std::cout << "step: " << config.step << '\n';
    std::cout << "population size: " << config.population_size << '\n';
    std::cout << "mutation rate: " << config.mutation_rate << '\n';
    std::cout << "cycles: " << config.cycles << '\n';
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


std::unique_ptr <std::vector <double>> generate_solution(const Func & function, const int multiplier, const int constant) {
    auto solution = std::make_unique <std::vector <double>> ();

    for (double x = 0; x <= 10.0; x += 0.1) {
        if (function == square) {
            solution->push_back(static_cast <double> (constant) + function(x) * static_cast <double> (multiplier));
        }
        else {
            solution->push_back(static_cast <double> (constant) + function(x * static_cast <double> (multiplier)));
        }
    }

    return solution;
}

std::unique_ptr <Individual> generate_solution(std::mt19937 & generator, const double min_value, const double max_value) {
    std::uniform_real_distribution <> distribution(min_value, max_value);

    auto solution = std::make_unique <std::vector <double>> (101);

    for (auto & element : * solution) {
        element = distribution(generator);
    }

    return std::make_unique <Individual> (0, std::move(solution));
}

void test_generated_solutions() {
    Config config = parse_config("../config.ini");

    std::random_device random_device;
    std::mt19937 generator(random_device());

    auto function = generate_solution(string_to_function_map[config.function], config.argument_multiplier, config.constant);

    for (const auto & a : * function) {
        std::cout << a << ' ';
    }
    std::cout << '\n';

    auto individual = generate_solution(generator, 0.0, 1.0);

    for (const auto & a : * (individual->solution)) {
        std::cout << a << ' ';
    }
    std::cout << '\n';
}

std::unique_ptr <std::vector <std::unique_ptr <Individual>>> create_population(const Config & config, std::mt19937 & generator) {
    auto population = std::make_unique <std::vector <std::unique_ptr <Individual>>> ();

    for (int i = 0; i < config.population_size; ++ i) {
        population->push_back(generate_solution(generator, -10.0, 10.0));
    }

    return population;
}

void ga_loop(const std::unique_ptr <std::vector <std::unique_ptr <Individual>>> & population, const Config & config, std::mt19937 & generator) {
    for (int cycle = 0; cycle < config.cycles; ++ cycle) {
        // step 1 crossover

        // step 2 mutation

        // step 3 fitness value evaluation

        // step 4 selection

        // break conditions

    }

    for (const auto & individual : * population) {
        for (int i = 0; i < 3; ++ i) {
            std::cout << individual->solution->at(i) << ' ';
        }

        std::cout << '\n';
    }

    std::cout << "all " << config.cycles << " cycles completed terminating the program...\n";
}


int main() {

    Config config = parse_config("../config.ini");
    print_config(config);

    std::random_device random_device;
    std::mt19937 generator(random_device());

    auto population = create_population(config, generator);

    ga_loop(population, config, generator);

    return 0;

}