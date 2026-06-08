#include <volt/centrosymmetry_engine.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <spdlog/spdlog.h>

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <mwm_csp.h>

namespace Volt{

using namespace Volt::Particles;

CentroSymmetryEngine::CentroSymmetryEngine(
    ParticleProperty* positions,
    const SimulationCell& cell,
    int numNeighbors,
    CentroSymmetryAnalysis::CSPMode mode
) : _positions(positions),
    _cell(cell),
    _k(numNeighbors),
    _mode(mode),
    _numBins(100),
    _histBinSize(1.0),
    _maxCSP(0.0)
{
    if(!_positions) throw std::runtime_error("CSP Engine: positions is null");
    if(_k < 2) throw std::runtime_error("CSP Engine: numNeighbors must be >= 2");
    if(_k > CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS) throw std::runtime_error("CSP Engine: numNeighbors too large");
    if((_k % 2) != 0) throw std::runtime_error("CSP Engine: numNeighbors must be even");
}

void CentroSymmetryEngine::perform(){
    const size_t n = _positions->size();

    if(n == 0){
        _csp.reset();
        _histCounts.reset();
        _maxCSP = 0.0;
        _histBinSize = 1.0;
        return;
    }

    // Output CSP per particle
    _csp = std::make_shared<ParticleProperty>(n, DataType::Double, 1, 0, true);

    // Build a kd-tree over the atoms once (numNeighbors == _k so each query
    // returns the _k nearest). Replaces the previous O(N^2) brute-force search.
    NearestNeighborFinder finder(_k);
    if(!finder.prepare(_positions, _cell)){
        for(size_t i = 0; i < n; i++) _csp->setDouble(i, 0.0);
        _maxCSP = 0.0;
        buildHistogram();
        return;
    }
    _finder = &finder;

    // Per-particle CSP in parallel. Each atom writes _csp[i] at a distinct
    // index and findKNearest builds its own thread-local kd-tree query, so
    // there is no shared mutable state and no data race.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
        [this](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i) computeParticleCSP(i);
        });

    _maxCSP = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n), 0.0,
        [this](const tbb::blocked_range<size_t>& r, double m){
            for(size_t i = r.begin(); i < r.end(); ++i) m = std::max(m, _csp->getDouble(i));
            return m;
        },
        [](double a, double b){ return std::max(a, b); });

    _finder = nullptr;
    buildHistogram();
}

void CentroSymmetryEngine::computeParticleCSP(size_t i){
    std::vector<Neighbor> neigh;
    neigh.reserve((size_t) _k);
    
    findKNearest(i, neigh);

    if((int) neigh.size() < _k){
        _csp->setDouble(i, 0.0);
        return;
    }

    double csp = computeCSPFromNeighbors(neigh);
    if(!std::isfinite(csp) || csp < 0) csp = 0.0;

    _csp->setDouble(i, csp);
}

void CentroSymmetryEngine::findKNearest(size_t i, std::vector<Neighbor>& out) const{
    out.clear();

    // kd-tree query: O(log N) per atom instead of O(N). The finder keeps the
    // _k nearest (constructed with numNeighbors == _k) and results() are sorted
    // ascending by distance. delta = neighbor - center, the same convention as
    // the previous minimum-image computation.
    NearestNeighborFinder::Query<CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS> q(*_finder);
    q.findNeighbors(i);

    const auto& res = q.results();
    out.reserve((size_t) res.size());
    for(int j = 0; j < res.size(); ++j){
        out.push_back(Neighbor{ res[j].distanceSq, res[j].delta });
    }
}

double CentroSymmetryEngine::computeCSPFromNeighbors(const std::vector<Neighbor>& neigh) const{
    const int numNN = (int) neigh.size();

    if(_mode == CentroSymmetryAnalysis::ConventionalMode){
        std::vector<double> pairs;
        pairs.reserve((size_t) numNN * (numNN - 1) / 2);

        for(int a = 0; a < numNN; a++){
            for(int b = a + 1; b < numNN; b++){
                Vector3 s = neigh[a].delta + neigh[b].delta;
                pairs.push_back((double) s.squaredLength());
            }
        }

        const int m = numNN / 2;
        if((int) pairs.size() < m) return 0.0;

        std::partial_sort(pairs.begin(), pairs.begin() + m, pairs.end());

        double csp = std::accumulate(pairs.begin(), pairs.begin() + m, 0.0);
        return csp;
    }else{
        // MatchingMode (MWM CSP)
        static_assert(CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS <= MWM_CSP_MAX_POINTS, "CSP neighbor limit mismatch");
        double P[CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS][3];
        for(int i = 0; i < numNN; i++){
            P[i][0] = (double)neigh[i].delta.x();
            P[i][1] = (double)neigh[i].delta.y();
            P[i][2] = (double)neigh[i].delta.z();
        }

        return (double)calculate_mwm_csp(numNN, P);
    }
}

void CentroSymmetryEngine::buildHistogram(){
    _histBinSize = (_maxCSP > 0.0) ? (1.01 * _maxCSP / (double)_numBins) : 0.0;
    if(_histBinSize <= 0.0) _histBinSize = 1.0;

    _histCounts = std::make_shared<ParticleProperty>(_numBins, DataType::Int64, 1, 0, true);
    for(size_t i = 0; i < _numBins; i++){
        _histCounts->setInt64(i, 0);
    }

    const size_t n = _positions->size();
    for(size_t i = 0; i < n; i++){
        const double v = _csp->getDouble(i);
        int bin = (int)(v / _histBinSize);
        if(bin < 0) bin = 0;
        if(bin >= (int)_numBins) continue;
        _histCounts->setInt64((size_t)bin, _histCounts->getInt64((size_t)bin) + 1);
    }
}


}
