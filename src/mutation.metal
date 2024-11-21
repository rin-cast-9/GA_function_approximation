#include <metal_stdlib>
using namespace metal;

kernel void mutate_kernel(
    device float * mutation_candidates [[ buffer(0) ]],
    device const float * mutation_values [[ buffer(1) ]],
    uint id [[ thread_position_in_grid ]]
) {
    mutation_candidates[id] += mutation_values[id];
}
