#include "common.h"
#include <cmath>
#include <unordered_map>
#include <vector>

#define BIN_SIZE 0.01

// Grid storage
std::vector<std::vector<particle_t*>> grid;
int num_bins;  // Number of bins in one dimension

// Function to compute 1D index from (x, y) bin coordinates
int bin_index(int x, int y) {
    return y * num_bins + x;
}

// Assign particles to bins
void bin_particles(particle_t* parts, int num_parts, double size) {
    num_bins = (int)(size / BIN_SIZE);  
    grid.clear();
    grid.resize(num_bins * num_bins);  // Allocate bins

    for (int i = 0; i < num_parts; ++i) {
        int x = (int)(parts[i].x / BIN_SIZE);
        int y = (int)(parts[i].y / BIN_SIZE);
        grid[bin_index(x, y)].push_back(&parts[i]);
    }
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

// Compute forces using binning
void compute_forces(particle_t* parts, int num_parts, double size) {
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;

        // Compute bin location
        int x = (int)(parts[i].x / BIN_SIZE);
        int y = (int)(parts[i].y / BIN_SIZE);

        // Iterate over neighboring bins
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int nx = x + dx;
                int ny = y + dy;

                // Ensure bin is within bounds
                if (nx >= 0 && nx < num_bins && ny >= 0 && ny < num_bins) {
                    for (particle_t* other : grid[bin_index(nx, ny)]) {
                        apply_force(parts[i], *other);
                    }
                }
            }
        }
    }
}

// Move particles
void move_particles(particle_t* parts, int num_parts, double size) {
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
}

// Simulation step
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    bin_particles(parts, num_parts, size);  // Assign to bins
    compute_forces(parts, num_parts, size); // Compute interactions
    move_particles(parts, num_parts, size); // Move particles
}

