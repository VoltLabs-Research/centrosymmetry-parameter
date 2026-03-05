#include <volt/centrosymmetry_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
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

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
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

    auto csp = engine.cspProperty();
    auto hist = engine.histogramCounts();

    // Build histogram array
    json histArray = json::array();
    if(hist){
        for(std::size_t b = 0; b < engine.numHistogramBins(); b++){
            histArray.push_back(hist->getInt64(b));
        }
    }

    json result;
    result["main_listing"] = {
        { "total_atoms", frame.natoms },
        { "histogram_bins", engine.numHistogramBins() },
        { "histogram_bin_size", engine.histogramBinSize() },
        { "max_csp", engine.maxCSP() },
        { "histogram_counts", histArray },
        { "histogram_interval", {
            { "start", 0.0 },
            { "end", engine.histogramBinSize() * (double)engine.numHistogramBins() }
        }}
    };

    json perAtom = json::array();
    for(int i = 0; i < frame.natoms; i++){
        perAtom.push_back({
            { "id", frame.ids[i] },
            { "csp", csp ? csp->getDouble(i) : 0.0 }
        });
    }
    result["per-atom-properties"] = perAtom;

    if(!outputBase.empty()){
        const std::string outputPath = outputBase + "_centrosymmetry.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Centrosymmetry msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write centrosymmetry msgpack: {}", outputPath);
        }
    }

    return result;
}

}

