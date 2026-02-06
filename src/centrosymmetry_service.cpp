#include <volt/centrosymmetry_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <spdlog/spdlog.h>

namespace Volt {

using namespace Volt::Particles;

CentroSymmetryService::CentroSymmetryService()
    : _k(12),
      _mode(CentroSymmetryAnalysis::ConventionalMode) {}

void CentroSymmetryService::setNumNeighbors(int k){
    _k = k;
}

void CentroSymmetryService::setMode(CentroSymmetryAnalysis::CSPMode mode){
    _mode = mode;
}

json CentroSymmetryService::compute(const LammpsParser::Frame& frame, const std::string& outputBase){
    auto start = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    if(_k < 2)
        return AnalysisResult::failure("numNeighbors must be >= 2");
    if(_k > CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS)
        return AnalysisResult::failure("numNeighbors too large");
    if((_k % 2) != 0)
        return AnalysisResult::failure("numNeighbors must be even");

    auto positions = FrameAdapter::createPositionProperty(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    spdlog::info("Computing CSP: k={}, mode={}", _k, (int)_mode);

    CentroSymmetryEngine engine(
        positions.get(),
        frame.simulationCell,
        _k,
        _mode
    );

    engine.perform();

    json result = AnalysisResult::success();
    AnalysisResult::addTiming(result, start);
    result["num_neighbors"] = _k;
    result["mode"] = (_mode == CentroSymmetryAnalysis::ConventionalMode) ? "conventional" : "matching";
    result["histogram_bins"] = engine.numHistogramBins();
    result["histogram_bin_size"] = engine.histogramBinSize();
    result["max_csp"] = engine.maxCSP();

    if(!outputBase.empty()){
        // TODO: Implement msgpack export in standalone package
        spdlog::warn("File output not yet implemented in standalone package");
    }

    return result;
}

}
