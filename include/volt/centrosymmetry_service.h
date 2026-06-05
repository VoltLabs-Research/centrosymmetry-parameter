#pragma once

#include <volt/centrosymmetry_engine.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt {

using json = nlohmann::json;

class CentrosymmetryService {
public:
    void setNumNeighbors(int k);
    void setMode(std::string mode);

    json compute(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    int _k = 12;
    CentroSymmetryAnalysis::CSPMode _mode = CentroSymmetryAnalysis::ConventionalMode;
};

} // namespace Volt
