#include <volt/cli/common.h>
#include <volt/centrosymmetry_service.h>
#include <oneapi/tbb/global_control.h>

using namespace Volt;
using namespace Volt::CLI;

static void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Centrosymmetry parameter (CSP)");
    std::cerr
        << "  --numNeighbors <int>          Even integer, <= 32. [default: 12]\n"
        << "  --mode <conventional|matching> [default: conventional]\n"
        << "  --threads <int>               Max worker threads. [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]){
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    const int requestedThreads = std::max(1, getInt(opts, "--threads", std::thread::hardware_concurrency() > 0
        ? static_cast<int>(std::thread::hardware_concurrency())
        : 1));
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(requestedThreads)
    );
    initLogging("volt-centrosymmetry");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    int k = getInt(opts, "--numNeighbors", 12);

    std::string modeStr = getString(opts, "--mode", "conventional");
    CentroSymmetryAnalysis::CSPMode mode = CentroSymmetryAnalysis::ConventionalMode;
    if(modeStr == "conventional") mode = CentroSymmetryAnalysis::ConventionalMode;
    else if(modeStr == "matching") mode = CentroSymmetryAnalysis::MatchingMode;
    else {
        spdlog::error("Invalid --mode. Use conventional or matching.");
        return 1;
    }

    CentroSymmetryService analyzer;
    analyzer.setNumNeighbors(k);
    analyzer.setMode(mode);

    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)){
        spdlog::error("CSP failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("CSP done. max_csp={}, bin_size={}",
        result.value("max_csp", 0.0),
        result.value("histogram_bin_size", 1.0));

    return 0;
}
