#pragma once
#include <cstdint>
#include <vector>

// Type definitions for consistent use across the library
typedef uint32_t index_t;   // For vertex and edge indices
typedef uint16_t comp_t;    // For component indices (using uint16_t to match your code examples)

typedef float real_t;       // For numerical computations
typedef std::vector<real_t> Vec3d;  // 3D point representation