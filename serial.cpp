#include "common.h"
#include <cmath>
#include <unordered_map>
#include <vector>

static double bin_size = cutoff; // bin size equals cutoff
static int num_bins_x, num_bins_y;
static std::vector<std::vector<std::vector<int>>> bins;

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

void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_size = cutoff;
    num_bins_x = static_cast<int>(size / bin_size) + 1;
    num_bins_y = static_cast<int>(size / bin_size) + 1;
    bins.resize(num_bins_x, std::vector<std::vector<int>>(num_bins_y));
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Clear the bins
    for (auto& row : bins) {
        for (auto& bin : row) {
            bin.clear();
        }
    }

    // Assign the particles to bins
    for (int i = 0; i < num_parts; ++i) {
        int bin_x = static_cast<int>(parts[i].x / bin_size);
        int bin_y = static_cast<int>(parts[i].y / bin_size);
        bin_x = std::max(0, std::min(bin_x, num_bins_x - 1));
        bin_y = std::max(0, std::min(bin_y, num_bins_y - 1));
        bins[bin_x][bin_y].push_back(i);
    }

    // Reset acceleration
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = 0.0;
        parts[i].ay = 0.0;
    }

    // Calculate forces
    for (int x = 0; x < num_bins_x; ++x) {
        for (int y = 0; y < num_bins_y; ++y) {
            auto& current_bin = bins[x][y];
            size_t size_current = current_bin.size();

            // Interactions within the current bin
            for (size_t i = 0; i < size_current; ++i) {
                for (size_t j = i + 1; j < size_current; ++j) {
                    int pi = current_bin[i];
                    int pj = current_bin[j];
                    apply_force(parts[pi], parts[pj]);
                    apply_force(parts[pj], parts[pi]);
                }
            }

            // Interactions with right bin
            if (x + 1 < num_bins_x) {
                auto& right_bin = bins[x + 1][y];
                for (int pi : current_bin) {
                    for (int pj : right_bin) {
                        apply_force(parts[pi], parts[pj]);
                        apply_force(parts[pj], parts[pi]);
                    }
                }
            }

            // Interactions with bottom bin
            if (y + 1 < num_bins_y) {
                auto& bottom_bin = bins[x][y + 1];
                for (int pi : current_bin) {
                    for (int pj : bottom_bin) {
                        apply_force(parts[pi], parts[pj]);
                        apply_force(parts[pj], parts[pi]);
                    }
                }
            }

            // Interactions with bottom-right bin
            if (x + 1 < num_bins_x && y + 1 < num_bins_y) {
                auto& br_bin = bins[x + 1][y + 1];
                for (int pi : current_bin) {
                    for (int pj : br_bin) {
                        apply_force(parts[pi], parts[pj]);
                        apply_force(parts[pj], parts[pi]);
                    }
                }
            }
            
            // Interactions with bottom-left bin
            if (x - 1 >= 0 && y + 1 < num_bins_y) {
                auto& bl_bin = bins[x - 1][y + 1];
                for (int pi : current_bin) {
                    for (int pj : bl_bin) {
                        apply_force(parts[pi], parts[pj]);
                        apply_force(parts[pj], parts[pi]);
                    }
                }
            }
        }
    }

    // Move particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}
