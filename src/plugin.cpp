#include <volt/plugin/plugin_entry.h>
#include <volt/centrosymmetry_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = CentrosymmetryService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--num_neighbors", "Even integer, <= 32", 12, &S::setNumNeighbors),
    opt("--mode", "conventional|matching", "conventional", &S::setMode),
};

VOLT_SERVICE_PLUGIN("volt-centrosymmetry", "Centrosymmetry Parameter (CSP)", S, bindings)
