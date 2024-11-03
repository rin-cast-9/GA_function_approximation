#include "GeneticAlgorithm.hpp"
#include "Individual.hpp"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <random>
#include <utility>
#include <vector>
#include <iomanip>
#include <iostream>
#include <future>


CrossoverFunction crossover_functions[] = {
    single_thread_crossover,
    parallel_crossover
};

MutationFunction mutation_functions[] = {
    execute_mutation,
    parallel_mutation
};

FitnessFunction fitness_functions[] = {
    evaluate_fitness,
    parallel_evaluate_fitness
};

std::unique_ptr <std::vector <double>> generate_solution(
    const Func & function,
    const int multiplier,
    const int constant,
    const double step
) {
    auto solution = std::make_unique <std::vector <double>> ();

    for (double x = 0; x <= 10.0; x += step) {
        if (function == square) {
            solution->push_back(static_cast <double> (constant) + function(x) * static_cast <double> (multiplier));
        }
        else {
            solution->push_back(static_cast <double> (constant) + function(x * static_cast <double> (multiplier)));
        }
    }

    return solution;
}

std::unique_ptr <Individual> generate_solution(
    std::mt19937 & generator,
    const double min_value,
    const double max_value,
    const double step
) {
    std::uniform_real_distribution <> distribution(min_value, max_value);

    auto solution = std::make_unique <std::vector <double>> ((int)(10.0 / step) + 1);

    for (auto & element : * solution) {
        element = distribution(generator);
    }

    return std::make_unique <Individual> (0, std::move(solution));
}

void crossover_type1(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
) {
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

void crossover_type2(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
) {
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

void crossover_type3(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <double> & child1_solution,
    std::vector <double> & child2_solution,
    std::mt19937 & generator
) {
    std::uniform_real_distribution <double> random_threshold(0.0, 1.0);
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

void crossover_type4(
    const Individual & parent1,
    const Individual & parent2,
    std::vector<double> & child1_solution,
    std::vector<double> & child2_solution,
    std::mt19937 & generator
) {
    size_t size = parent1.solution->size();

    auto & parent1_solution = * parent1.solution;
    auto & parent2_solution = * parent2.solution;

    for (size_t i = 0; i < size; ++ i) {
        std::uniform_real_distribution <double> distribution(
            std::min(parent1_solution[i], parent2_solution[i]),
            std::max(parent1_solution[i], parent2_solution[i])
        );

        child1_solution[i] = distribution(generator);
        child2_solution[i] = distribution(generator);
    }
}

std::pair<std::unique_ptr <Individual>, std::unique_ptr <Individual>> execute_crossover(
    const Individual & parent1,
    const Individual & parent2,
    std::mt19937 & generator,
    const int crossover_strategy
) {
    const size_t solution_size = parent1.solution->size();
    auto child_individual1 = std::make_unique <std::vector <double>> (solution_size);
    auto child_individual2 = std::make_unique <std::vector <double>> (solution_size);

    switch (crossover_strategy) {
    case 0:
        crossover_type1(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    case 1:
        crossover_type2(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    case 2:
        crossover_type3(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    case 3:
        crossover_type4(parent1, parent2, * child_individual1, * child_individual2, generator);
        break;
    }

    auto child1 = std::make_unique <Individual> (0.0, std::move(child_individual1));
    auto child2 = std::make_unique <Individual> (0.0, std::move(child_individual2));

    return std::make_pair(std::move(child1), std::move(child2));
}

void parallel_crossover(
    std::vector <std::unique_ptr <Individual>> & population,
    std::mt19937 & generator,
    const int population_size,
    const int max_population_size,
    const int crossover_strategy
) {
    // std::vector <std::pair <std::unique_ptr <Individual>, std::unique_ptr <Individual>>> new_individuals;

    population.resize(max_population_size);

    size_t pairs_to_crossover = (max_population_size - population_size) / 2;

    std::vector <std::future <std::pair <std::unique_ptr <Individual>, std::unique_ptr <Individual>>>> futures;

    for (size_t i = 0; i < pairs_to_crossover; ++ i) {
        auto crossover_candidates = choose_crossover_candidates(generator, population_size);

        futures.push_back(std::async(std::launch::async, [&](){
            return execute_crossover(
                * population[crossover_candidates.first],
                * population[crossover_candidates.second],
                generator,
                crossover_strategy
            );
        }));
    }

    for (size_t i = 0; i < futures.size(); ++ i) {
        auto children = futures[i].get();

        population[population_size + (i * 2)] = std::move(children.first);
        population[population_size + (i * 2) + 1] = std::move(children.second);
    }
}

void single_thread_crossover(
    std::vector <std::unique_ptr<Individual>> & population,
    std::mt19937 & generator,
    const int population_size,
    const int max_population_size,
    const int crossover_strategy
) {
    while (population.size() < max_population_size) {
        auto crossover_candidates = choose_crossover_candidates(generator, population_size);

        auto children = execute_crossover(
            * population[crossover_candidates.first],
            * population[crossover_candidates.second],
            generator,
            crossover_strategy
        );

        population.push_back(std::move(children.first));
        population.push_back(std::move(children.second));
    }
}

std::pair <int, int> choose_crossover_candidates(
    std::mt19937 & generator,
    const int population_size
) {
    std::geometric_distribution <int> distribution(0.2);

    int candidate1 = distribution(generator) % population_size;
    int candidate2 = distribution(generator) % population_size;

    do {
        candidate2 = distribution(generator) % population_size;
    } while (candidate1 == candidate2);

    return {candidate1, candidate2};
}

std::vector <size_t> choose_mutation_candidate(
    std::mt19937 & generator,
    const size_t max_population_size,
    const double mutation_rate
) {
    auto candidates = std::vector <size_t> ();

    std::bernoulli_distribution distribution(mutation_rate);

    for (size_t i = 0; i < max_population_size; ++ i) {
        if (distribution(generator)) {
            candidates.push_back(i);
        }
    }

    return candidates;
}

void execute_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const double average_fitness_value,
    const double min_range,
    const double max_range,
    const double step,
    std::mt19937 & generator
) {
    const double range = std::abs(max_range - min_range);
    const double domain = (((int)(10.0 / step) + 1) * range);
    const double mutation_strength = (average_fitness_value < (domain * 0.03)) ? (average_fitness_value / domain / 10) : (average_fitness_value / domain / 3);

    std::uniform_real_distribution <double> distribution(-range, range);

    for (const auto & candidate_index : mutation_candidates) {
        auto & individual = population[candidate_index];
        auto & solution = * individual->solution;

        for (double & element : solution) {
            double mutation_value = distribution(generator) * mutation_strength;
            element += mutation_value;
        }
    }
}

void parallel_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const double average_fitness_value,
    const double min_range,
    const double max_range,
    const double step,
    std::mt19937 & generator
) {
    const double range = std::abs(max_range - min_range);
    const double mutation_strength = average_fitness_value / (((int)(10.0 / step) + 1) * range) / 10.0;

    std::uniform_real_distribution <double> distribution(-range, range);

    std::vector <std::future <void>> futures;

    for (const auto & candidate_index : mutation_candidates) {
        futures.push_back(std::async(std::launch::async, [&, candidate_index]() {
            auto & individual = population[candidate_index];
            auto & solution = * individual->solution;

            for (double & element : solution) {
                double mutation_value = distribution(generator) * mutation_strength;
                element += mutation_value;
            }
        }));
    }

    for (auto & future : futures) {
        future.get();
    }
}

double evaluate_average_fitness(
    const std::vector <std::unique_ptr <Individual>> & population
) {
    return std::accumulate(
        population.begin(),
        population.end(),
        0.0,
        [](double sum, const std::unique_ptr<Individual> & individual) {
            return sum + individual->fitness;
        }
    ) / (double) population.size();
}

double evaluate_maximum_fitness_value(
    const std::vector <std::unique_ptr <Individual>> & population
) {
    return std::min_element(
        population.begin(),
        population.end(),
        [](const std::unique_ptr <Individual> & a, const std::unique_ptr <Individual> & b) {
            return a->fitness < b->fitness;
        }
    )->get()->fitness;
}

void print_cycle_data(
    const int cycle,
    const double average_fitness_value,
    const double max_fitness_value,
    const int mutation_amount
) {
    std::cout << "\033[31;1m# " << std::setw(7) << cycle << "\033[0m \033[32;1mAVG: " << std::setw(8) << average_fitness_value << "\033[0m \033[34;1mMAX: " << std::setw(8) << max_fitness_value << "\033[0m \033[35;1mMUT: " << std::setw(3) << mutation_amount << "\033[0m\n";
}

void print_graph_data_to_file(
    const std::string & file_path,
    const std::vector <double> & target,
    const std::vector <double> & best_solution,
    const double step
) {
    std::ofstream out_file(file_path);

    if (!out_file) {
        return;
    }

    out_file << step << '\n';
    for (const auto & e : target) {
        out_file << e << ' ';
    }
    out_file << '\n';
    for (const auto & e : best_solution) {
        out_file << e << ' ';
    }
    out_file << '\n';
}

void ga_loop(
    const std::unique_ptr <std::vector <std::unique_ptr <Individual>>> & population,
    const std::unique_ptr <std::vector <double>> & target,
    const Config & config,
    std::mt19937 & generator
) {
    fitness_functions[config.computation_mode](
        * target,
        * population
    );

    double average_fitness_value = evaluate_average_fitness(* population);
    double max_fitness_value = evaluate_maximum_fitness_value(* population);

    for (int cycle = 0; cycle < config.cycles; ++ cycle) {

        perform_selection(* population, config.population_size);

        // step 1 crossover
        if (population->size() < 2) {
            std::cerr << "population size is too small for crossover.\n";
            throw std::runtime_error("population size must be at least 2 for crossover.\n");
        }

        crossover_functions[config.computation_mode](
            * population,
            generator,
            config.population_size,
            config.max_population_size,
            config.crossover_strategy
        );

        // print_population(* population);

        // step 2 mutation
        auto mutation_candidates = choose_mutation_candidate(generator, config.max_population_size, config.mutation_rate);

        mutation_functions[config.computation_mode](
            * population,
            mutation_candidates,
            average_fitness_value,
            config.min_range,
            config.max_range,
            config.step,
            generator
        );

        // step 3 fitness value evaluation
        fitness_functions[config.computation_mode](
            * target,
            * population
        );

        average_fitness_value = evaluate_average_fitness(* population);
        max_fitness_value = evaluate_maximum_fitness_value(* population);

        if (cycle % 100 == 0) {
            print_cycle_data(cycle, average_fitness_value, max_fitness_value, mutation_candidates.size());
            // print_population(* population);
        }

        // step 4 selection
        // perform_selection(* population, config.population_size);

        // break conditions
        if (max_fitness_value < 0.5) {
            break;
        }
    }

    print_graph_data_to_file(std::string(PROJECT_ROOT_DIR) + "/graph_data.txt", * target, * (population->at(0)->solution), config.step);

    std::cout << "\nall " << config.cycles << " cycles completed terminating the program...\n";
}
