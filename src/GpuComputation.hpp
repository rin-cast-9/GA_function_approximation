#pragma once

#ifndef GPU_COMPUTATION_HPP
#define GPU_COMPUTATION_HPP

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#include "Metal/MTLComputePipeline.hpp"
#include "Metal/MTLLibrary.hpp"
#include "Individual.hpp"
#include <memory>
#include <vector>

namespace MTL {
    class Device;
    class CommandQueue;
    class Buffer;
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
);

void gpu_evaluate_fitness(
    const std::vector<float> & target,
    const std::vector<std::unique_ptr<Individual>> & population,
    const std::vector<float> & flattened_population,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
);

void gpu_mutate(
    const std::vector<std::unique_ptr<Individual>> & population,
    const std::vector<float> & flattened_solutions,
    const std::vector<float> & mutation_values,
    const std::vector<size_t> & mutation_candidates,
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library,
    MTL::ComputePipelineState * & pipeline_state
);

void initialize_metal_resources(
    MTL::Device * & device,
    MTL::CommandQueue * & command_queue,
    MTL::Library * & library,
    const std::string & metallib_path
);

void release_metal_resources(
    MTL::Device * device,
    MTL::CommandQueue * command_queue,
    MTL::Library * library
);

void release_pipeline_state(
    MTL::ComputePipelineState * pipeline_state
);

#endif
