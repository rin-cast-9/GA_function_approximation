#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include "GpuComputation.hpp"
#include "Individual.hpp"
#include <cstddef>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include <stdexcept>

void initialize_metal_resources(
    MTL::Device * & device,
    MTL::CommandQueue * & command_queue,
    MTL::Library * & library,
    const std::string & metallib_path
) {
    device = MTL::CreateSystemDefaultDevice();

    if (!device) {
        std::cerr << "Metal device not found or not supported.\n";
        throw std::runtime_error("Metal device not found or not supported.\n");
    }

    command_queue = device->newCommandQueue();

    if (!command_queue) {
        std::cerr << "Failed to create command queue.\n";
        throw std::runtime_error("Failed to create command queue.\n");
    }

    NS::Error * error = nullptr;
    NS::String * library_path = NS::String::string((std::string(PROJECT_BINARY_DIR) + "/default.metallib").c_str(), NS::UTF8StringEncoding);


    library = device->newLibrary(library_path, & error);
    if (!library || error) {
        std::cerr << "Failed to load library: " << error->localizedDescription()->utf8String() << '\n';
        throw std::runtime_error("Failed to load library.");
    }
}

void release_metal_resources(
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library
) {
    command_queue->release();
    library->release();
    device->release();
}

void release_pipeline_state(
    MTL::ComputePipelineState * pipeline_state
) {
    pipeline_state->release();
}

void gpu_crossover(
    const std::pair <std::unique_ptr<std::vector<float>>, std::unique_ptr<std::vector<float>>> & flatten_solutions,
    std::vector <std::unique_ptr<Individual>> & population,
    const int population_size,
    const int max_population_size,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
) {
    const size_t individual_size = population[0]->solution->size();
    const int pairs = (max_population_size - population_size) / 2;

//    std::vector <float> child1_data(pairs * individual_size);
//    std::vector <float> child2_data(pairs * individual_size);

    MTL::CommandBuffer * command_buffer = command_queue->commandBuffer();
    MTL::ComputeCommandEncoder * encoder = command_buffer->computeCommandEncoder();

    size_t total_elements = pairs * individual_size;

    MTL::Buffer * buffer_parent1 = device->newBuffer(flatten_solutions.first->data(), sizeof(float) * total_elements, MTL::ResourceStorageModeShared);
    MTL::Buffer * buffer_parent2 = device->newBuffer(flatten_solutions.second->data(), sizeof(float) * total_elements, MTL::ResourceStorageModeShared);
    MTL::Buffer * buffer_child1 = device->newBuffer(sizeof(float) * total_elements, MTL::ResourceStorageModeShared);
    MTL::Buffer * buffer_child2 = device->newBuffer(sizeof(float) * total_elements, MTL::ResourceStorageModeShared);

    std::vector <int> random_selector(total_elements);
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution <float> distribution(0.0, 1.0);

    for (auto & value : random_selector) {
        value = distribution(generator) < 0.5 ? 0 : 1;
    }
    MTL::Buffer * buffer_random_selector = device->newBuffer(random_selector.data(), sizeof(int) * random_selector.size(), MTL::ResourceStorageModeShared);

    if (!pipeline_state) {
        NS::Error * error = nullptr;
        MTL::Function * function = library->newFunction(NS::String::string("crossover_type3_kernel", NS::UTF8StringEncoding));
        if (!function) {
            std::cerr << "Failed to load kernel function: crossover_type3_kernel.\n";
            throw std::runtime_error("Failed to load kernel function: crossover_type3_kernel.\n");
        }

        pipeline_state = device->newComputePipelineState(function, & error);
    }

    encoder->setComputePipelineState(pipeline_state);
    encoder->setBuffer(buffer_parent1, 0, 0);
    encoder->setBuffer(buffer_parent2, 0, 1);
    encoder->setBuffer(buffer_child1, 0, 2);
    encoder->setBuffer(buffer_child2, 0, 3);
    encoder->setBuffer(buffer_random_selector, 0, 4);

    MTL::Size grid_size = MTL::Size((max_population_size - population_size) / 2 * individual_size, 1, 1);
    MTL::Size thread_group_size = MTL::Size(1, 1, 1);
    encoder->dispatchThreads(grid_size, thread_group_size);

    encoder->endEncoding();
    command_buffer->commit();
    command_buffer->waitUntilCompleted();

    float * child1_raw = static_cast<float *>(buffer_child1->contents());
    float * child2_raw = static_cast<float *>(buffer_child2->contents());

    for (size_t i = 0; i < pairs; ++ i) {
        auto child1 = std::make_unique<std::vector<float>>(child1_raw + i * individual_size, child1_raw + (i + 1) * individual_size);

        auto child2 = std::make_unique<std::vector<float>>(child2_raw + i * individual_size, child2_raw + (i + 1) * individual_size);

        population.push_back(std::make_unique<Individual>(0.0, std::move(child1)));
        population.push_back(std::make_unique<Individual>(0.0, std::move(child2)));
    }

    buffer_parent1->release();
    buffer_parent2->release();
    buffer_child1->release();
    buffer_child2->release();
    buffer_random_selector->release();
    encoder->release();
    command_buffer->release();

}

void gpu_evaluate_fitness(
    const std::vector<float> & target,
    const std::vector<std::unique_ptr<Individual>> & population,
    const std::vector<float> & flattened_population,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
) {
    const size_t individual_size = target.size();
    const size_t target_padding_needed = (4 - individual_size % 4) % 4;

    std::vector<float> input_data; // target + population
    input_data.reserve(individual_size + target_padding_needed + flattened_population.size());
    input_data.insert(input_data.end(), target.begin(), target.end());
    input_data.insert(input_data.end(), target_padding_needed, 0.0f);
    input_data.insert(input_data.end(), flattened_population.begin(), flattened_population.end());

    MTL::CommandBuffer * command_buffer = command_queue->commandBuffer();
    MTL::ComputeCommandEncoder * encoder = command_buffer->computeCommandEncoder();

    MTL::Buffer * input_data_buffer = device->newBuffer(input_data.data(), sizeof(float) * input_data.size(), MTL::ResourceStorageModeShared);
    MTL::Buffer * fitness_buffer = device->newBuffer(sizeof(float) * input_data.size(), MTL::ResourceStorageModeShared);
    MTL::Buffer * output_sum_buffer = device->newBuffer(sizeof(float) * population.size(), MTL::ResourceStorageModeShared);

    if (!pipeline_state) {
        NS::Error * error = nullptr;
        MTL::Function * function = library->newFunction(NS::String::string("evaluate_fitness_kernel", NS::UTF8StringEncoding));
        if (!function) {
            std::cerr << "Failed to load kernel function: evaluate_fitness_kernel.\n";
            throw std::runtime_error("Failed to load kernel function: evaluate_fitness_kernel.\n");
        }

        pipeline_state = device->newComputePipelineState(function, & error);
    }

    const size_t total_padded_size = individual_size + target_padding_needed;

    encoder->setComputePipelineState(pipeline_state);
    encoder->setBuffer(input_data_buffer, 0, 0);
    encoder->setBuffer(fitness_buffer, 0, 1);
    encoder->setBytes(& total_padded_size, sizeof(total_padded_size), 2);
    encoder->setBuffer(output_sum_buffer, 0, 3);

    MTL::Size grid_size = MTL::Size(total_padded_size, population.size() + 1, 1);
    MTL::Size thread_group_size = MTL::Size(4, 1, 1);
    size_t ldata_size = sizeof(float) * thread_group_size.width;
    encoder->setThreadgroupMemoryLength(ldata_size, 0);
    encoder->dispatchThreads(grid_size, thread_group_size);

    encoder->endEncoding();
    command_buffer->commit();
    command_buffer->waitUntilCompleted();

    float * output_sum_data = static_cast<float *>(output_sum_buffer->contents());
    for (size_t i = 0; i < population.size(); ++ i) {
        population[i]->fitness = output_sum_data[i];
    }

    input_data_buffer->release();
    fitness_buffer->release();
    output_sum_buffer->release();
    encoder->release();
    command_buffer->release();
}


void gpu_mutate(
    const std::vector<std::unique_ptr<Individual>> & population,
    const std::vector<float> & flattened_solutions,
    const std::vector<float> & mutation_values,
    const std::vector<size_t> & mutation_candidates,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
) {
    MTL::CommandBuffer * command_buffer = command_queue->commandBuffer();
    MTL::ComputeCommandEncoder * encoder = command_buffer->computeCommandEncoder();

    MTL::Buffer * mutation_candidates_buffer = device->newBuffer(flattened_solutions.data(), sizeof(float) * flattened_solutions.size(), MTL::ResourceStorageModeShared);
    MTL::Buffer * mutation_values_buffer = device->newBuffer(mutation_values.data(), sizeof(float) * mutation_values.size(), MTL::ResourceStorageModeShared);

    if (!pipeline_state) {
        NS::Error * error = nullptr;
        MTL::Function * function = library->newFunction(NS::String::string("mutate_kernel", NS::UTF8StringEncoding));
        if (!function) {
            std::cerr << "Failed to load kernel function: mutate_kernel.\n";
            throw std::runtime_error("Failed to load kernel function: mutate_kernel.\n");
        }

        pipeline_state = device->newComputePipelineState(function, & error);
        std::cout << pipeline_state->maxTotalThreadsPerThreadgroup() << '\n';
    }

    encoder->setComputePipelineState(pipeline_state);
    encoder->setBuffer(mutation_candidates_buffer, 0, 0);
    encoder->setBuffer(mutation_values_buffer, 0, 1);

    MTL::Size grid_size = MTL::Size(flattened_solutions.size(), 1, 1);
    MTL::Size thread_group_size = MTL::Size(1, 1, 1);
    encoder->dispatchThreads(grid_size, thread_group_size);

    encoder->endEncoding();
    command_buffer->commit();
    command_buffer->waitUntilCompleted();

    float * mutated_candidates = static_cast<float *>(mutation_candidates_buffer->contents());
    const size_t individual_size = population[0]->solution->size();

    for (size_t i = 0; i < mutation_candidates.size(); ++ i) {
        const size_t index_in_population = mutation_candidates[i];

        auto mutated_individual = std::make_unique<std::vector<float>>(mutated_candidates + i * individual_size, mutated_candidates + (i + 1) * individual_size);

        population[index_in_population]->solution = std::move(mutated_individual);
    }

    mutation_candidates_buffer->release();
    mutation_values_buffer->release();
    encoder->release();
    command_buffer->release();
}
