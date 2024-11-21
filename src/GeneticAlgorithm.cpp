#include "GeneticAlgorithm.hpp"
#include "GpuComputation.hpp"
#include "Individual.hpp"
#include "Metal/MTLComputePipeline.hpp"
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
    parallel_crossover,
    // gpu_crossover_setup
};

MutationFunction mutation_functions[] = {
    execute_mutation,
    parallel_mutation,
    parallel_mutation
};

FitnessFunction fitness_functions[] = {
    evaluate_fitness,
    parallel_evaluate_fitness,
    // parallel_evaluate_fitness
};

std::unique_ptr <std::vector <float>> generate_solution(
    const Func & function,
    const int multiplier,
    const int constant,
    const double step
) {
    auto solution = std::make_unique <std::vector <float>> ();

    for (double x = 0; x <= 10.0; x += step) {
        if (function == square) {
            solution->push_back(static_cast <float> (constant) + function(x) * static_cast <float> (multiplier));
        }
        else {
            solution->push_back(static_cast <float> (constant) + function(x * static_cast <float> (multiplier)));
        }
    }

    return solution;
}

std::unique_ptr <Individual> generate_solution(
    std::mt19937 & generator,
    const float min_value,
    const float max_value,
    const double step
) {
    std::uniform_real_distribution <> distribution(min_value, max_value);

    auto solution = std::make_unique <std::vector <float>> ((int)(10.0 / step) + 1);

    for (auto & element : * solution) {
        element = distribution(generator);
    }

    return std::make_unique <Individual> (0, std::move(solution));
}

void crossover_type1(
    const Individual & parent1,
    const Individual & parent2,
    std::vector <float> & child1_solution,
    std::vector <float> & child2_solution,
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
    std::vector <float> & child1_solution,
    std::vector <float> & child2_solution,
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
    std::vector <float> & child1_solution,
    std::vector <float> & child2_solution,
    std::mt19937 & generator
) {
    std::uniform_real_distribution <float> random_threshold(0.0, 1.0);
    float threshold = random_threshold(generator);
    size_t size = parent1.solution->size();

    for (size_t i = 0; i < size; ++ i) {
        float y = random_threshold(generator);
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
    std::vector<float> & child1_solution,
    std::vector<float> & child2_solution,
    std::mt19937 & generator
) {
    size_t size = parent1.solution->size();

    auto & parent1_solution = * parent1.solution;
    auto & parent2_solution = * parent2.solution;

    for (size_t i = 0; i < size; ++ i) {
        std::uniform_real_distribution <float> distribution(
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
    auto child_individual1 = std::make_unique <std::vector <float>> (solution_size);
    auto child_individual2 = std::make_unique <std::vector <float>> (solution_size);

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

std::pair <std::unique_ptr <std::vector <float>>, std::unique_ptr <std::vector <float>>> flatten_solutions(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <std::pair<int, int>> & chosen_crossover_candidates
) {
    auto flattened_solutions_parent1 = std::make_unique <std::vector <float>> ();
    auto flattened_solutions_parent2 = std::make_unique <std::vector <float>> ();
    flattened_solutions_parent1->reserve(chosen_crossover_candidates.size());
    flattened_solutions_parent2->reserve(chosen_crossover_candidates.size());

    for (const auto & candidate_pair : chosen_crossover_candidates) {
        const auto & parent1_solution = * population[candidate_pair.first]->solution;
        const auto & parent2_solution = * population[candidate_pair.second]->solution;

        flattened_solutions_parent1->insert(flattened_solutions_parent1->end(), parent1_solution.begin(), parent1_solution.end());
        flattened_solutions_parent2->insert(flattened_solutions_parent2->end(), parent2_solution.begin(), parent2_solution.end());
    }

    return std::make_pair(std::move(flattened_solutions_parent1), std::move(flattened_solutions_parent2));
}

std::unique_ptr<std::vector<float>> flatten_population(
    const std::vector<std::unique_ptr<Individual>> & population
) {
    auto flattened_population = std::make_unique <std::vector<float>> ();

    size_t solution_size = population[0]->solution->size();
    size_t padding_needed = (4 - solution_size % 4) % 4;

    flattened_population->reserve(population.size() * (solution_size + padding_needed));

    for (const auto & individual : population) {
        flattened_population->insert(flattened_population->end(),
                                     individual->solution->begin(),
                                     individual->solution->end());
        flattened_population->insert(flattened_population->end(), padding_needed, 0.0f);
    }

    return flattened_population;
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

void gpu_crossover_setup(
    std::vector <std::unique_ptr <Individual>> & population,
    std::mt19937 & generator,
    const int population_size,
    const int max_population_size,
    const int crossover_strategy,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
) {
    std::vector<std::pair<int, int>> crossover_candidates_batch;

    const int pairs = (max_population_size - population_size) / 2;

    for (int i = 0; i < pairs; ++ i) {
        auto crossover_candidates = choose_crossover_candidates(generator, population_size);
        crossover_candidates_batch.push_back(crossover_candidates);
    }

    auto flattened_solutions = flatten_solutions(population, crossover_candidates_batch);

    gpu_crossover(flattened_solutions, population, population_size, max_population_size, device, command_queue, library, pipeline_state);
}

void gpu_fitness_setup(
    const std::vector <float> & target,
    const std::vector <std::unique_ptr<Individual>> & population,
    const int population_size,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
) {
    auto flattened_population = flatten_population(population);

    gpu_evaluate_fitness(target, population, * flattened_population, device, command_queue, library, pipeline_state);
}

void gpu_mutation_setup(const std::vector<std::unique_ptr<Individual>> & population,
                        const std::vector<size_t> & mutation_candidates,
                        const float average_fitness_value,
                        const float min_range,
                        const float max_range,
                        const double step,
                        std::mt19937 & generator,
                        MTL::Device * device,
                        MTL::CommandQueue * command_queue,
                        MTL::Library * library,
                        MTL::ComputePipelineState * & pipeline_state) {
    auto flattened_solutions = flatten_solutions(population, mutation_candidates);

    auto mutation_values = generate_mutation_values(flattened_solutions->size(), average_fitness_value, min_range, max_range, step, generator);

    gpu_mutate(population, * flattened_solutions, * mutation_values, mutation_candidates, device, command_queue, library, pipeline_state);
}

std::unique_ptr<std::vector<float>> flatten_solutions(const std::vector<std::unique_ptr<Individual>> & population,
                                                      const std::vector<size_t> & mutation_candidates) {
    const size_t individual_size = population[0]->solution->size();
    const size_t mutation_candidates_size = mutation_candidates.size();

    auto flattened_mutation_candidates = std::make_unique<std::vector<float>>();
    flattened_mutation_candidates->reserve(individual_size * mutation_candidates_size);

    for (const auto & mutation_candidate_index : mutation_candidates) {
        auto mutation_candidate = * population[mutation_candidate_index]->solution;

        flattened_mutation_candidates->insert(flattened_mutation_candidates->end(), mutation_candidate.begin(), mutation_candidate.end());
    }

    return flattened_mutation_candidates;
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

std::unique_ptr<std::vector <size_t>> choose_mutation_candidate(
    std::mt19937 & generator,
    const size_t max_population_size,
    const float mutation_rate
) {
    auto candidates = std::make_unique<std::vector <size_t>>();

    std::bernoulli_distribution distribution(mutation_rate);

    for (size_t i = 0; i < max_population_size; ++ i) {
        if (distribution(generator)) {
            candidates->push_back(i);
        }
    }

    return candidates;
}

std::unique_ptr<std::vector<float>> generate_mutation_values(const size_t mutation_values_needed,
                                                             const float average_fitness_value,
                                                             const float min_range,
                                                             const float max_range,
                                                             const double step,
                                                             std::mt19937 & generator) {
    auto mutation_values = std::make_unique<std::vector<float>>();
    mutation_values->resize(mutation_values_needed);

    const float range = std::abs(max_range - min_range);
    const float domain = (((int)(10.0 / step) + 1) * range);
    const float mutation_strength = (average_fitness_value / domain) * 0.5;

    std::uniform_real_distribution<float> distribution(-range, range);

    for (auto & mutation_value : * mutation_values) {
        mutation_value = distribution(generator) * mutation_strength;
    }

    return mutation_values;
}

void execute_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const float average_fitness_value,
    const float min_range,
    const float max_range,
    const double step,
    std::mt19937 & generator
) {
    const float range = std::abs(max_range - min_range);
    const float domain = (((int)(10.0 / step) + 1) * range);
    // const float mutation_strength = (average_fitness_value < (domain * 0.03)) ? (average_fitness_value / domain / 10) : (average_fitness_value / domain / 3);
    const float mutation_strength = (average_fitness_value / domain) * 0.5;

    std::uniform_real_distribution <float> distribution(-range, range);

    for (const auto & candidate_index : mutation_candidates) {
        auto & individual = population[candidate_index];
        auto & solution = * individual->solution;

        for (float & element : solution) {
            float mutation_value = distribution(generator) * mutation_strength;
            element += mutation_value;
        }
    }
}

void parallel_mutation(
    const std::vector <std::unique_ptr <Individual>> & population,
    const std::vector <size_t> & mutation_candidates,
    const float average_fitness_value,
    const float min_range,
    const float max_range,
    const double step,
    std::mt19937 & generator
) {
    const float range = std::abs(max_range - min_range);
    const float mutation_strength = average_fitness_value / (((int)(10.0 / step) + 1) * range) / 10.0;

    std::uniform_real_distribution <float> distribution(-range, range);

    std::vector <std::future <void>> futures;

    for (const auto & candidate_index : mutation_candidates) {
        futures.push_back(std::async(std::launch::async, [&, candidate_index]() {
            auto & individual = population[candidate_index];
            auto & solution = * individual->solution;

            for (float & element : solution) {
                float mutation_value = distribution(generator) * mutation_strength;
                element += mutation_value;
            }
        }));
    }

    for (auto & future : futures) {
        future.get();
    }
}

float evaluate_average_fitness(
    const std::vector <std::unique_ptr <Individual>> & population
) {
    return std::accumulate(
        population.begin(),
        population.end(),
        0.0,
        [](float sum, const std::unique_ptr<Individual> & individual) {
            return sum + individual->fitness;
        }
    ) / (float) population.size();
}

float evaluate_maximum_fitness_value(
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
    const float average_fitness_value,
    const float max_fitness_value,
    const int mutation_amount
) {
   std::cout << "\033[31;1m# " << std::setw(7) << cycle << "\033[0m \033[32;1mAVG: " << std::setw(8) << average_fitness_value << "\033[0m \033[34;1mMAX: " << std::setw(8) << max_fitness_value << "\033[0m \033[35;1mMUT: " << std::setw(3) << mutation_amount << "\033[0m\n";
}

void print_graph_data_to_file(
    const std::string & file_path,
    const std::vector <float> & target,
    const std::vector <float> & best_solution,
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

    std::cout << target.size() << '\n';
    std::cout << best_solution.size() << '\n';
}

void ga_loop(
    const std::unique_ptr <std::vector <std::unique_ptr <Individual>>> & population,
    const std::unique_ptr <std::vector <float>> & target,
    const Config & config,
    std::mt19937 & generator
) {

    MTL::Device * device = nullptr;
    MTL::CommandQueue * command_queue = nullptr;
    MTL::Library * library = nullptr;
    MTL::ComputePipelineState * pipeline_state_crossover = nullptr;
    MTL::ComputePipelineState * pipeline_state_fitness = nullptr;
    MTL::ComputePipelineState * pipeline_state_mutation = nullptr;

    const std::string metallib_path = "/default.metallib";

    if (config.computation_mode == 2) {
        initialize_metal_resources(device, command_queue, library, metallib_path);

        gpu_fitness_setup(* target, * population, config.population_size, device, command_queue, library, pipeline_state_fitness);
    }
    else {
        fitness_functions[config.computation_mode](
            * target,
            * population
        );
    }

    float average_fitness_value = evaluate_average_fitness(* population);
    float max_fitness_value = evaluate_maximum_fitness_value(* population);

    for (int cycle = 0; cycle < config.cycles; ++ cycle) {

        perform_selection(* population, config.population_size);

        // step 1 crossover
        if (population->size() < 2) {
            std::cerr << "population size is too small for crossover.\n";
            throw std::runtime_error("population size must be at least 2 for crossover.\n");
        }

        if (config.computation_mode == 2) {
            gpu_crossover_setup(* population, generator, config.population_size, config.max_population_size, config.crossover_strategy, device, command_queue, library, pipeline_state_crossover);
        }
        else {
            crossover_functions[config.computation_mode](
                * population,
                generator,
                config.population_size,
                config.max_population_size,
                config.crossover_strategy
            );
        }

        // print_population(* population);

        // step 2 mutation
        auto mutation_candidates = choose_mutation_candidate(generator, config.max_population_size, config.mutation_rate);

        if (config.computation_mode == 2) {
            gpu_mutation_setup(* population, * mutation_candidates, average_fitness_value, config.min_range, config.max_range, config.step, generator, device, command_queue, library, pipeline_state_mutation);
        }
        else {
            mutation_functions[config.computation_mode](
                * population,
                * mutation_candidates,
                average_fitness_value,
                config.min_range,
                config.max_range,
                config.step,
                generator
            );
        }

        // step 3 fitness value evaluation
        if (config.computation_mode == 2) {
            gpu_fitness_setup(* target, * population, config.population_size, device, command_queue, library, pipeline_state_fitness);
        }
        else {
            fitness_functions[config.computation_mode](
                * target,
                * population
            );
        }

        average_fitness_value = evaluate_average_fitness(* population);
        max_fitness_value = evaluate_maximum_fitness_value(* population);

        if (cycle % config.print_interval == 0) {
            print_cycle_data(cycle, average_fitness_value, max_fitness_value, mutation_candidates->size());
            // print_population(* population);
        }

        // step 4 selection
        perform_selection(* population, config.population_size);

        // break conditions
        if (max_fitness_value < config.approximation_tolerance) {
            break;
        }
    }

    release_metal_resources(device, command_queue, library);
    release_pipeline_state(pipeline_state_crossover);
    release_pipeline_state(pipeline_state_fitness);
    release_pipeline_state(pipeline_state_mutation);

    print_graph_data_to_file(std::string(PROJECT_ROOT_DIR) + "/graph_data.txt", * target, * (population->at(0)->solution), config.step);
    std::cout << "data written!\n";
}
