#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// ============================================================================
// Decomposition : fonctions de scalarisation pour MOEA/D
// ============================================================================

// Classe abstraite
class Decomposition {
  public:
    virtual ~Decomposition() = default;

    // Fonction virtuelle pure
    virtual double value(const std::vector<double>& obj,
                          const std::vector<double>& weight,
                          const std::vector<double>& z,
                          const std::vector<double>& zmax) const = 0;
};


// Tchebycheff normalisé
class TchebycheffNormalized : public Decomposition {
  public:
    double value(const std::vector<double>& obj,
                            const std::vector<double>& weight,
                            const std::vector<double>& z,
                            const std::vector<double>& zmax) const override {
        double g = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < obj.size(); ++k) {
            double range = zmax[k] - z[k];
            if (range < 1e-12) range = 1e-12; // Évite la division par 0

            double diff = std::fabs(obj[k] - z[k]) / range * weight[k];
            g = std::max(g, diff);
        }
        return g;
    }
};


// PBI (Penalty-based Boundary Intersection)
class PBI : public Decomposition {
  public:
    explicit PBI(double theta = 5.0) : theta_(theta) {}

    double value(const std::vector<double>& obj,
                 const std::vector<double>& weight,
                 const std::vector<double>& z,
                 const std::vector<double>& zmax) const override {
        size_t M = obj.size();

        double normW = 0.0;
        for (size_t k = 0; k < M; ++k) normW += weight[k] * weight[k];
        normW = std::sqrt(normW);
        if (normW < 1e-12) normW = 1e-12;

        // d1
        double dot = 0.0;
        for (size_t k = 0; k < M; ++k) {
            double range = zmax[k] - z[k];
            if (range < 1e-12) range = 1e-12;

            double yk = (obj[k] - z[k]) / range; 
            dot += yk * weight[k];
        }
        double d1 = dot / normW;

        // d2
        double d2sq = 0.0;
        for (size_t k = 0; k < M; ++k) {
            double range = zmax[k] - z[k];
            if (range < 1e-12) range = 1e-12;

            double yk = (obj[k] - z[k]) / range;
            double proj = d1 * weight[k] / normW;
            double diff = yk - proj;
            d2sq += diff * diff;
        }
        double d2 = std::sqrt(d2sq);

        return d1 + theta_ * d2;
    }

  private:
    double theta_;
};
