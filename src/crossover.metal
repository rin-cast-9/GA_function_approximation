#include <metal_stdlib>

using namespace metal;

kernel void crossover_type3_kernel(
    device const float * parent1 [[ buffer(0) ]],
    device const float * parent2 [[ buffer(1) ]],
    device float * child1 [[ buffer(2) ]],
    device float * child2 [[ buffer(3) ]],
    device const int * random_values [[ buffer(4) ]],
    uint id [[ thread_position_in_grid ]]
) {
    child1[id] = random_values[id] * parent1[id] + (1 - (random_values[id])) * parent2[id];
    child2[id] = random_values[id] * parent2[id] + (1 - (random_values[id])) * parent1[id];
}
