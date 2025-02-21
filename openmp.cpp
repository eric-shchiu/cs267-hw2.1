#include "common.h"
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include <vector>

#define BIN_SIZE 0.01

static double bin_size = cutoff;
static int num_bins_x, num_bins_y;
static std::vector<std::vector<std::vector<int>>> bins;

// Modified apply_force with atomic operations
void apply_force(particle_t& particle, particle_t& neighbor) {
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef = (1 - cutoff / r) / r2 / mass;

    #pragma omp atomic
    particle.ax += coef * dx;
    #pragma omp atomic
    particle.ay += coef * dy;
}

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

// Process neighboring bins
void process_neighbor(const std::vector<int>& cur_bin, 
    const std::vector<int>& nb_bin,
    particle_t* parts) {
    for (int p1 : cur_bin) {
        for (int p2 : nb_bin) {
            apply_force(parts[p1], parts[p2]);
            apply_force(parts[p2], parts[p1]);
            }
    }
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_size = cutoff;
    num_bins_x = static_cast<int>(size / bin_size) + 1;
    num_bins_y = static_cast<int>(size / bin_size) + 1;
    bins.resize(num_bins_x, std::vector<std::vector<int>>(num_bins_y));
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Parallel bin clearing
    #pragma omp for collapse(2)
    for (int x = 0; x < num_bins_x; ++x) {
        for (int y = 0; y < num_bins_y; ++y) {
            bins[x][y].clear();
        }
    }

    // Single-threaded bin assignment
    #pragma omp single
    for (int i = 0; i < num_parts; ++i) {
        int bin_x = static_cast<int>(parts[i].x / bin_size);
        int bin_y = static_cast<int>(parts[i].y / bin_size);
        bin_x = std::max(0, std::min(bin_x, num_bins_x - 1));
        bin_y = std::max(0, std::min(bin_y, num_bins_y - 1));
        bins[bin_x][bin_y].push_back(i);
    }

    // Parallel acceleration reset
    #pragma omp for
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = 0.0;
        parts[i].ay = 0.0;
    }

    // Parallel force computation
    #pragma omp for collapse(2)
    for (int x = 0; x < num_bins_x; ++x) {
        for (int y = 0; y < num_bins_y; ++y) {
            auto& cur_bin = bins[x][y];
            const size_t n = cur_bin.size();

            // interact within the same bin
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    particle_t& p1 = parts[cur_bin[i]];
                    particle_t& p2 = parts[cur_bin[j]];
                    apply_force(p1, p2);
                    apply_force(p2, p1);
                }
            }

            // interact with right bin
            if (x+1 < num_bins_x) process_neighbor(cur_bin, bins[x+1][y], parts);
            
            // interact with bottom bin
            if (y+1 < num_bins_y) process_neighbor(cur_bin, bins[x][y+1], parts);
            
            // interact with bottom-right bin
            if (x+1 < num_bins_x && y+1 < num_bins_y) 
                process_neighbor(cur_bin, bins[x+1][y+1], parts);
            
            // interact with bottom-left bin
            if (x > 0 && y+1 < num_bins_y) 
                process_neighbor(cur_bin, bins[x-1][y+1], parts);
        }
    }

    // Parallel particle movement
    #pragma omp for
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}

