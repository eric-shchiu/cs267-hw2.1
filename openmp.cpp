#include "common.h"
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include <vector>

// Put any static global variables here that you will use throughout the simulation.

#define BIN_SIZE 0.00001  // Adjust based on interaction radius

// Hash function for cell indexing
struct CellHash {
    size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ std::hash<int>()(p.second);
    }
};

// Spatial grid: maps cell (x, y) to list of particles
std::unordered_map<std::pair<int, int>, std::vector<particle_t*>, CellHash> grid;

// Function to get cell index from position
std::pair<int, int> get_cell(double x, double y) {
    return { (int)(x / BIN_SIZE), (int)(y / BIN_SIZE) };
}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here
}

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // grid.clear();

    // Parallel Force Computation
    #pragma omp parallel for
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        std::pair<int, int> cell = get_cell(parts[i].x, parts[i].y);

        // Check forces from own and neighboring bins
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                std::pair<int, int> neighbor = { cell.first + dx, cell.second + dy };
                if (grid.find(neighbor) != grid.end()) {
                    for (particle_t* other : grid[neighbor]) {
                        apply_force(parts[i], *other);
                    }
                }
            }
        }
    }

    // Parallel Particle Movement
    #pragma omp parallel for
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}
