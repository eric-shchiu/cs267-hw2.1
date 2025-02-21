#include "common.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

// Memory Layout Optimization: SoA Structure Definition
struct Particles {
    std::vector<double> x, y;      // Location
    std::vector<double> vx, vy;    // Velocity
    std::vector<double> ax, ay;    // Acceleration

    explicit Particles(int num) {
        x.resize(num); y.resize(num);
        vx.resize(num); vy.resize(num);
        ax.resize(num); ay.resize(num);
    }
};

// Global split-box data structure (flattened)
static double bin_size = cutoff;      // bin size is fixed to cutoff
static int num_bins_x, num_bins_y;    // Number of sub-cases
static std::vector<int> bin_data;     // Stores particle indexes for all bins (one-dimensional contiguous memory)
static std::vector<int> bin_offsets;  // Starting offset of each bin
static std::vector<int> bin_counts;   // Number of particles per bin

// Force Calculation Functions (Adapted to SoA)
inline void apply_force_soa(
    double x1, double y1, double& ax1, double& ay1,
    double x2, double y2, double& ax2, double& ay2
) {
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    const double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff) return;

    const double r2_clamped = std::fmax(r2, min_r * min_r);
    const double r = std::sqrt(r2_clamped);
    const double coef = (1.0 - cutoff / r) / r2_clamped / mass;

    ax1 += coef * dx;
    ay1 += coef * dy;
    ax2 -= coef * dx;
    ay2 -= coef * dy;
}

// Move Functions (Adaptation to SoA)
inline void move_soa(
    double& x, double& y,
    double& vx, double& vy,
    double ax, double ay,
    double size
) {
    vx += ax * dt;
    vy += ay * dt;
    x += vx * dt;
    y += vy * dt;

    if (x < 0) {
        x = -x;
        vx = -vx;
    } else if (x > size) {
        x = 2 * size - x;
        vx = -vx;
    }

    if (y < 0) {
        y = -y;
        vy = -vy;
    } else if (y > size) {
        y = 2 * size - y;
        vy = -vy;
    }
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_size = cutoff;
    num_bins_x = static_cast<int>(size / bin_size) + 1;
    num_bins_y = static_cast<int>(size / bin_size) + 1;

    const int total_bins = num_bins_x * num_bins_y;
    bin_data.resize(num_parts);
    bin_offsets.resize(total_bins + 1);
    bin_counts.resize(total_bins, 0);
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    Particles soa_parts(num_parts);
    for (int i = 0; i < num_parts; ++i) {
        soa_parts.x[i] = parts[i].x;
        soa_parts.y[i] = parts[i].y;
        soa_parts.vx[i] = parts[i].vx;
        soa_parts.vy[i] = parts[i].vy;
        soa_parts.ax[i] = 0.0;
        soa_parts.ay[i] = 0.0;
    }

    // Empty the count
    std::fill(bin_counts.begin(), bin_counts.end(), 0);

    // Sub-case statistics
    for (int i = 0; i < num_parts; ++i) {
        int bin_x = static_cast<int>(soa_parts.x[i] / bin_size);
        bin_x = std::max(0, std::min(bin_x, num_bins_x - 1));
        int bin_y = static_cast<int>(soa_parts.y[i] / bin_size);
        bin_y = std::max(0, std::min(bin_y, num_bins_y - 1));
        const int bin_idx = bin_x * num_bins_y + bin_y;
        bin_counts[bin_idx]++;
    }

    // Compute the prefix sum
    bin_offsets[0] = 0;
    std::partial_sum(bin_counts.begin(), bin_counts.end(), bin_offsets.begin() + 1);

    // Filling data
    std::fill(bin_counts.begin(), bin_counts.end(), 0);
    for (int i = 0; i < num_parts; ++i) {
        int bin_x = static_cast<int>(soa_parts.x[i] / bin_size);
        bin_x = std::max(0, std::min(bin_x, num_bins_x - 1));
        int bin_y = static_cast<int>(soa_parts.y[i] / bin_size);
        bin_y = std::max(0, std::min(bin_y, num_bins_y - 1));
        const int bin_idx = bin_x * num_bins_y + bin_y;
        const int pos = bin_offsets[bin_idx] + bin_counts[bin_idx]++;
        bin_data[pos] = i;
    }

    // Calculate forces
    for (int i = 0; i < num_parts; ++i) {
        int bin_x = static_cast<int>(soa_parts.x[i] / bin_size);
        bin_x = std::max(0, std::min(bin_x, num_bins_x - 1));
        int bin_y = static_cast<int>(soa_parts.y[i] / bin_size);
        bin_y = std::max(0, std::min(bin_y, num_bins_y - 1));

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                const int nb_x = bin_x + dx;
                const int nb_y = bin_y + dy;
                if (nb_x < 0 || nb_x >= num_bins_x || nb_y < 0 || nb_y >= num_bins_y) continue;

                const int bin_idx = nb_x * num_bins_y + nb_y;
                const int start = bin_offsets[bin_idx];
                const int end = bin_offsets[bin_idx + 1];

                for (int k = start; k < end; ++k) {
                    const int j = bin_data[k];
                    if (j > i) {
                        apply_force_soa(
                            soa_parts.x[i], soa_parts.y[i], soa_parts.ax[i], soa_parts.ay[i],
                            soa_parts.x[j], soa_parts.y[j], soa_parts.ax[j], soa_parts.ay[j]
                        );
                    }
                }
            }
        }
    }

    // Move particles
    for (int i = 0; i < num_parts; ++i) {
        move_soa(
            soa_parts.x[i], soa_parts.y[i],
            soa_parts.vx[i], soa_parts.vy[i],
            soa_parts.ax[i], soa_parts.ay[i],
            size
        );
    }

    // Copy back
    for (int i = 0; i < num_parts; ++i) {
        parts[i].x = soa_parts.x[i];
        parts[i].y = soa_parts.y[i];
        parts[i].vx = soa_parts.vx[i];
        parts[i].vy = soa_parts.vy[i];
    }
}
