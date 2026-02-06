#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/centrosymmetry_engine.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt{
using json = nlohmann::json;

class CentroSymmetryService{
public:
    CentroSymmetryService();

    void setNumNeighbors(int k);
    void setMode(CentroSymmetryAnalysis::CSPMode mode);

    json compute(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    int _k;
    CentroSymmetryAnalysis::CSPMode _mode;
};

}
