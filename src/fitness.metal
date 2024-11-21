#include <metal_stdlib>

using namespace metal;

kernel void evaluate_fitness_kernel(
    device const float * data [[ buffer(0) ]],
    device float * fitness [[ buffer(1) ]],
    constant uint & individual_size [[ buffer(2) ]],
    device atomic_float * output_sum [[ buffer(3) ]],
    threadgroup float * ldata [[ threadgroup(0) ]],
    uint2 gid [[ thread_position_in_grid ]],
    uint2 lid [[ thread_position_in_threadgroup ]],
    uint2 lsize [[ threads_per_threadgroup ]]
) {
    uint row = gid.y;
    uint col = gid.x;
    uint index = row * individual_size + col;

    if (col > individual_size) {
        return;
    }

    fitness[index] = (row == 0) ? data[index] : fabs(data[col] - data[index]);

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float val = (row > 0 && col < individual_size) ? fitness[index] : 0.0f;

    val = quad_sum(val);

    if (quad_is_first()) {
        uint quadgroup_id = (lid.y * lsize.x + lid.x) / 4;
        ldata[quadgroup_id] = val;
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (lid.y == 0 && lid.x < lsize.x / 4) {
        val = ldata[lid.x];
        val = quad_sum(val);

        if (quad_is_first() && lid.x == 0 && row > 0) {
            fitness[row * individual_size] = val;
            atomic_fetch_add_explicit(& output_sum[row - 1], val, memory_order_relaxed);
        }
    }
}
