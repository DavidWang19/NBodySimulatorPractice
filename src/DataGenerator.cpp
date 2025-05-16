#include "SystemInitializers.h"
#include "SystemOutput.h"

int main() {
    std::vector<nbs::PhysicalObject> objects;
    nbs::Inhomogeneous2D::initialize(objects, 100000, 10, 200);
    nbs::SystemOutput output("../../../../initial/I2D_100k_200.txt");
    output.output_system_arrangement_raw(objects);
    return 0;
}