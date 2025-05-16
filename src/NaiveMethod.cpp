#include "EvolutionSimulation.h"
#include "NaiveMethod.h"
#include "SystemInitializers.h"
#include "SystemOutput.h"

#include <iostream>
#include <memory>
#include <vector>

using std::cout;
using std::endl;

int main() {
    std::vector<nbs::PhysicalObject> objects;
    nbs::ReadFromFile::initialize(objects, "../../../../initial/U3D_err.txt");
    std::unique_ptr<nbs::EvolutionAlgorithmInterface> method =
        std::make_unique<nbs::NaiveMethod>("../../../../results/U3D_err_naive.log");
    nbs::EvolutionSimulation simulation(std::move(objects), std::move(method), 2);
    nbs::SystemOutput output_util("../../../../results/U3D_err_naive.txt");
    output_util.output_system_arrangement(simulation);
    cout << "Outputting initial system arrangement" << endl;
    for (int i = 0; i < 3000; ++i) {
        simulation.evolve();
        cout << "Evolved system #" << i << endl;
        output_util.output_system_arrangement(simulation);
    }
    return 0;
}