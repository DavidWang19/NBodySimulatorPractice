#include "BarnesHutMethod.h"
#include "EvolutionSimulation.h"
#include "SystemInitializers.h"
#include "SystemOutput.h"

#include <iostream>
#include <memory>
#include <vector>

using std::cout;
using std::endl;

int main() {
    std::vector<nbs::PhysicalObject> objects;
    nbs::ReadFromFile::initialize(objects, "../../../../initial/I2D_100k_200.txt");
    std::unique_ptr<nbs::EvolutionAlgorithmInterface> method =
        std::make_unique<nbs::BarnesHutMethod>(200, 200, 1, 0.3,
                                               "../../../../results/I2D_100k_200_BH.log");
    static_cast<nbs::BarnesHutMethod*>(method.get())
        ->set_output_file_name("../../../../results/I2D_100k_200_octree.txt");
    nbs::EvolutionSimulation simulation(std::move(objects), std::move(method), 0.2);
    nbs::SystemOutput output_util("../../../../results/I2D_100k_200_BH.txt");
    output_util.output_system_arrangement(simulation);
    cout << "Outputting initial system arrangement" << endl;
    for (int i = 0; i < 3000; ++i) {
        simulation.evolve();
        cout << "Evolved system #" << i << endl;
        output_util.output_system_arrangement(simulation);
    }
    return 0;
}