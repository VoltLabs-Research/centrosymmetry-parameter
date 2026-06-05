#include <volt/plugin/plugin_main.h>
#include <volt/plugin/analysis_plugin.h>
#include <volt/plugin/output_serializer.h>
#include <volt/centrosymmetry_engine.h>

namespace Volt {

class CentrosymmetryPlugin : public Plugin::AnalysisPlugin<CentrosymmetryPlugin> {
public:
    int k = 12;
    CentroSymmetryAnalysis::CSPMode mode = CentroSymmetryAnalysis::ConventionalMode;

    std::string validate(const LammpsParser::Frame&) {
        if (k < 2) return "numNeighbors must be >= 2";
        if (k > CentroSymmetryAnalysis::MAX_CSP_NEIGHBORS) return "numNeighbors too large";
        if ((k % 2) != 0) return "numNeighbors must be even";
        return "";
    }

    json run(const LammpsParser::Frame& frame,
             const std::shared_ptr<Particles::ParticleProperty>& positions,
             const std::string&) {
        CentroSymmetryEngine engine(positions.get(), frame.simulationCell, k, mode);
        engine.perform();

        _csp = engine.cspProperty();
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
        return result;
    }

    void serialize(const LammpsParser::Frame& frame,
                   const std::shared_ptr<Particles::ParticleProperty>&,
                   const json& result,
                   const std::string& outputBase) {
        auto csp = _csp;
        Plugin::serializePluginOutput(outputBase, frame, result, {
            .summaryFileSuffix = "_centrosymmetry",
            .bucketResolver = [](std::size_t) { return std::string("All"); },
            .atomFieldWriter = [csp](MsgpackWriter& w, std::size_t i, int& count) {
                count = 1;
                w.write_key("csp");
                w.write_double(csp ? csp->getDouble(i) : 0.0);
            }
        });
    }

private:
    std::shared_ptr<Particles::ParticleProperty> _csp;
};

} // namespace Volt

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-centrosymmetry",
    .description = "Centrosymmetry parameter (CSP)",
    .options = {
        {"--num_neighbors", "int", "Even integer, <= 32", "12"},
        {"--mode", "string", "conventional|matching", "conventional"},
    }
};

VOLT_PLUGIN_MAIN(descriptor, [](const auto& opts, const Volt::LammpsParser::Frame& frame,
                                 const Volt::LammpsParser::Frame*, const std::string& outputBase) {
    Volt::CentrosymmetryPlugin plugin;
    plugin.k = Volt::CLI::getInt(opts, "--num_neighbors", 12);
    std::string modeStr = Volt::CLI::getString(opts, "--mode", "conventional");
    if (modeStr == "matching")
        plugin.mode = Volt::CentroSymmetryAnalysis::MatchingMode;
    else if (modeStr != "conventional") {
        spdlog::error("Invalid --mode. Use conventional or matching.");
        return Volt::AnalysisResult::failure("Invalid --mode value");
    }
    return plugin.compute(frame, outputBase);
})
