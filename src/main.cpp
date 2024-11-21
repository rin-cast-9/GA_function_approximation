#include "Individual.hpp"
#include "Config.hpp"
#include "GeneticAlgorithm.hpp"
#include "Utilities.hpp"

#include <iostream>
#include <string>
#include <map>
#include <random>
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

void test_generated_solutions() {
    Config config = parse_config("../config.conf");

    std::random_device random_device;
    std::mt19937 generator(random_device());

    auto function = generate_solution(string_to_function_map[config.function], config.argument_multiplier, config.constant, config.step);

    for (const auto & a : * function) {
        std::cout << a << ' ';
    }
    std::cout << '\n';

    auto individual = generate_solution(generator, 0.0, 1.0, config.step);

    for (const auto & a : * (individual->solution)) {
        std::cout << a << ' ';
    }
    std::cout << '\n';
}

int main() {

    Config config = parse_config(std::string(PROJECT_ROOT_DIR) + "/config.conf");
    print_config(config);

    std::random_device random_device;
    std::mt19937 generator(random_device());

    auto population = create_population(config, generator);
    auto target = generate_solution(string_to_function_map[config.function], config.argument_multiplier, config.constant, config.step);

    ga_loop(population, target, config, generator);

    return 0;

}
