#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


double square(const double x) {
    return std::pow(x, 2);
}

using Func = double (*)(double);

std::map <std::string, Func> StringToFunctionMap = {
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
    int constant;
    std::string function;
    int argument_multiplier;
};

Config parseConfig(const std::string & filePath) {
    Config config;
    boost::property_tree::ptree ptree;
    boost::property_tree::ini_parser::read_ini(filePath, ptree);

    config.step = ptree.get <double> ("step");
    config.population_size = ptree.get <int> ("population_size");
    config.mutation_rate = ptree.get <double> ("mutation_rate");
    config.constant = ptree.get <int> ("function.constant");
    config.function = ptree.get <std::string> ("function.function");
    config.argument_multiplier = ptree.get <int> ("function.argument_multiplier");

    return config;
}

void printConfig(const Config & config) {
    std::cout << "step: " << config.step << '\n';
    std::cout << "population size: " << config.population_size << '\n';
    std::cout << "mutation rate: " << config.mutation_rate << '\n';
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


int main() {

    std::cout << "hello, ga program kokokara desu!\n";

    Config config = parseConfig("../config.ini");
    printConfig(config);

    return 0;

}