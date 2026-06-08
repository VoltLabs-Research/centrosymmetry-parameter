#include <volt/centrosymmetry_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/plugin/output_serializer.h>
#include <spdlog/spdlog.h>

namespace Volt {

void CentrosymmetryService::setNumNeighbors(int k) { _k = k; }

void CentrosymmetryService::setMode(std::string mode) {
    if (mode == "matching")
        _mode = CentroSymmetryAnalysis::MatchingMode;
    else
        _mode = CentroSymmetryAnalysis::ConventionalMode;
}

json CentrosymmetryService::compute(const LammpsParser::Frame& frame, const std::string& outputBase) {
    if (frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");
    if (_k < 2) return AnalysisResult::failure("numNeighbors must be >= 2");
    if (_k > CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS) return AnalysisResult::failure("numNeighbors too large");
    if ((_k % 2) != 0) return AnalysisResult::failure("numNeighbors must be even");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if (!positions) return AnalysisResult::failure("Failed to create position property");

    CentroSymmetryEngine engine(positions.get(), frame.simulationCell, _k, _mode);
    engine.perform();

    auto csp = engine.cspProperty();
    auto hist = engine.histogramCounts();

    json histogramRows = json::array();
    if (hist) {
        for (std::size_t b = 0; b < engine.numHistogramBins(); ++b) {
            const double start = engine.histogramBinSize() * static_cast<double>(b);
            histogramRows.push_back({
                {"bin", static_cast<int>(b)},
                {"start", start},
                {"end", start + engine.histogramBinSize()},
                {"count", hist->getInt64(b)}
            });
        }
    }

    const double histogramEnd = engine.histogramBinSize() * static_cast<double>(engine.numHistogramBins());
    json result;
    result["main_listing"] = {
        {"total_atoms", frame.natoms},
        {"histogram_bins", engine.numHistogramBins()},
        {"histogram_bin_size", engine.histogramBinSize()},
        {"max_csp", engine.maxCSP()},
        {"histogram_start", 0.0},
        {"histogram_end", histogramEnd}
    };
    result["sub_listings"] = {{"histogram", histogramRows}};

    if (!outputBase.empty()) {
        Plugin::serializePluginOutput(outputBase, frame, result, {
            .summaryFileSuffix = "_centrosymmetry",
            .bucketResolver = [](std::size_t) { return std::string("All"); },
            .atomFieldWriter = [&csp](MsgpackWriter& w, std::size_t i, int& count) {
                count = 1;
                w.write_key("csp");
                w.write_double(csp ? csp->getDouble(i) : 0.0);
            },
            .perAtomFieldWriter = [&csp](MsgpackWriter& w, std::size_t i, int& count) {
                count = 1;
                w.write_key("csp");
                w.write_double(csp ? csp->getDouble(i) : 0.0);
            }
        });
    }

    return result;
}

} // namespace Volt
