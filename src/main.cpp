#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <memory>
#include <set>
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
    int max_population_size;
    double mutation_rate;
    int cycles;
    int crossover_strategy;
    int constant;
    std::string function;
    int argument_multiplier;
};

struct Individual {
    double fitness;
    std::unique_ptr <std::vector <double>> solution;

    Individual(double fit, std::unique_ptr <std::vector <double>> sol) : fitness(fit), solution(std::move(sol)) {}
};

void validate_config_step(double step) {
    std::set <double> valid_steps = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};

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

void validate_mutation_rate(double mutation_rate) {
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
    if (crossover_strategy < 0 || crossover_strategy > 2) {
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

Config parse_config(const std::string & file_path) {
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

    validate_config_step(config.step);
    validate_population_size(config.population_size, config.max_population_size);
    validate_mutation_rate(config.mutation_rate);
    validate_cycles(config.cycles);
    validate_crossover_strategy(config.crossover_strategy);
    validate_function(config.constant, config.function, config.argument_multiplier);

    return config;
}

void print_config(const Config & config) {
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
    Config config = parse_config("../config.conf");

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

void crossover_type1(const Individual & parent1, const Individual & parent2, std::vector <double> & child1_solution, std::vector <double> & child2_solution, std::mt19937 & generator) {
    std::uniform_int_distribution <int> distribution(0, parent1.solution->size() - 1);
    int crossover_point = distribution(generator);
    size_t size = parent1.solution->size();

    for (size_t i = 0; i < crossover_point; ++ i) {
        child1_solution[i] = parent1.solution->at(i);
        child2_solution[i] = parent2.solution->at(i);
    }
    for (size_t i = crossover_point; i < parent2.solution->size(); ++ i) {
        child1_solution[i] = parent2.solution->at(i);
        child2_solution[i] = parent1.solution->at(i);
    }
}

void crossover_type2(const Individual & parent1, const Individual & parent2, std::vector <double> & child1_solution, std::vector <double> & child2_solution, std::mt19937 & generator) {
    std::uniform_int_distribution <int> distribution(0, 1);
    size_t size = parent1.solution->size();

    for (size_t i = 0; i < size; ++ i) {
        if (distribution(generator) == 0) {
            child1_solution[i] = parent1.solution->at(i);
            child2_solution[i] = parent2.solution->at(i);
        }
        else {
            child1_solution[i] = parent2.solution->at(i);
            child2_solution[i] = parent1.solution->at(i);
        }
    }
}

void crossover_type3(const Individual & parent1, const Individual & parent2, std::vector <double> & child1_solution, std::vector <double> & child2_solution, std::mt19937 & generator) {
    std::uniform_int_distribution <double> random_threshold(0.0, 1.0);
    double threshold = random_threshold(generator);
    size_t size = parent1.solution->size();

    for (size_t i = 0; i < size; ++ i) {
        double y = random_threshold(generator);
        if (y < threshold) {
            child1_solution[i] = parent1.solution->at(i);
            child2_solution[i] = parent2.solution->at(i);
        }
        else {
            child1_solution[i] = parent2.solution->at(i);
            child2_solution[i] = parent1.solution->at(i);
        }
    }
}

std::pair<std::unique_ptr <Individual>, std::unique_ptr <Individual>> execute_crossover(const Individual & parent1, const Individual & parent2, std::mt19937 & generator, int crossover_strategy) {
    auto child_individual1 = std::make_unique <std::vector <double>> (parent1.solution->size());
    auto child_individual2 = std::make_unique <std::vector <double>> (parent1.solution->size());
    
    switch (crossover_strategy) {
    case 1:
        crossover_type1(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    case 2:
        crossover_type2(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    case 3:
        crossover_type3(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    }

    auto child1 = std::make_unique <Individual> (0.0, std::move(child_individual1));
    auto child2 = std::make_unique <Individual> (0.0, std::move(child_individual2));

    return std::make_pair(std::move(child1), std::move(child2));
}

std::pair <int, int> choose_crossover_candidates(std::mt19937 & generator, int population_size) {
    std::uniform_int_distribution <int> distribution(0, population_size - 1);

    int candidate1 = distribution(generator);
    int candidate2 = distribution(generator);

    do {
        candidate2 = distribution(generator);
    } while (candidate1 == candidate2);

    return {candidate1, candidate2};
}

void ga_loop(const std::unique_ptr <std::vector <std::unique_ptr <Individual>>> & population, const Config & config, std::mt19937 & generator) {
    for (int cycle = 0; cycle < config.cycles; ++ cycle) {
        // step 1 crossover
        if (population->size() < 2) {
            std::cerr << "population size is too small for crossover.\n";
            throw std::runtime_error("population size must be at least 2 for crossover.\n");
        }

        while (population->size() < config.max_population_size) {
            auto crossover_candidates = choose_crossover_candidates(generator, config.population_size);

            auto children = execute_crossover(
                * population->at(crossover_candidates.first),
                * population->at(crossover_candidates.second),
                generator,
                config.crossover_strategy
            );

            population->push_back(std::move(children.first));
            population->push_back(std::move(children.second));
        }

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

    Config config = parse_config("../config.conf");
    print_config(config);

    std::random_device random_device;
    std::mt19937 generator(random_device());

    // auto population = create_population(config, generator);

    // ga_loop(population, config, generator);

    return 0;

}