#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


struct Config {
    double step;
    std::string function;
    int population_size;
    double mutation_rate;
};

Config parseConfig(const std::string & filePath) {
    Config config;
    boost::property_tree::ptree ptree;
    boost::property_tree::ini_parser::read_ini(filePath, ptree);

    config.step = ptree.get <double> ("step");
    config.function = ptree.get <std::string> ("function");
    config.population_size = ptree.get <int> ("population_size");
    config.mutation_rate = ptree.get <double> ("mutation_rate");

    return config;
}

void printConfig(const Config & config) {
    std::cout << "step: " << config.step << '\n';
    std::cout << "function: " << config.function << '\n';
    std::cout << "population size: " << config.population_size << '\n';
    std::cout << "mutation rate: " << config.mutation_rate << '\n';
}


int main() {

    std::cout << "hello, ga program kokokara desu!\n";

    Config config = parseConfig("../config.ini");
    printConfig(config);

    return 0;

}