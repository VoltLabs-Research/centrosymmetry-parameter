#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <memory>
#include <vector>
#include <cstddef>

namespace Volt{

using namespace Particles;

class CentroSymmetryAnalysis{
public:
    enum CSPMode{
        ConventionalMode = 0,
        MatchingMode = 1
    };

    enum { MAX_CSP_NEIGHBORS = 32 };
};

class CentroSymmetryEngine{
public:
    CentroSymmetryEngine(
        ParticleProperty* positions,
        const SimulationCell& cell,
        int numNeighbors,
        CentroSymmetryAnalysis::CSPMode mode
    );

    void perform();

    std::shared_ptr<ParticleProperty> cspProperty() const{
        return _csp;
    }

    std::shared_ptr<ParticleProperty> histogramCounts() const{
        return _histCounts;
    }

    double histogramBinSize() const{
        return _histBinSize;
    }

    size_t numHistogramBins() const{
        return _numBins;
    }

    double maxCSP() const{
        return _maxCSP;
    }

private:
    struct Neighbor{
        double dist2;
        Vector3 delta;
    };

    void computeParticleCSP(size_t i);
    double computeCSPFromNeighbors(const std::vector<Neighbor>& neigh) const;

    void findKNearest(size_t i, std::vector<Neighbor>& out) const;

    void buildHistogram();

    ParticleProperty* _positions;
    SimulationCell _cell;
    int _k;
    CentroSymmetryAnalysis::CSPMode _mode;

    std::shared_ptr<ParticleProperty> _csp;
    std::shared_ptr<ParticleProperty> _histCounts;

    size_t _numBins;
    double _histBinSize;
    double _maxCSP;
};

}
