#include "SystemOutput.h"

using namespace nbs;

void SystemOutput::output_system_arrangement(EvolutionSimulation& simulation) {
    double time = simulation.get_time();
    m_output << "Time: " << time << "\n";
    m_output << "No.,mass,x,y,z,v_x,v_y,v_z\n";
    output_system_arrangement_raw(simulation.get_objects());
    m_output << "--------------------------------\n";
}

void SystemOutput::output_system_arrangement_raw(const std::vector<PhysicalObject>& objects) {
    for (size_t i = 0; i < objects.size(); ++i) {
        const auto& object = objects[i];
        const auto& position = object.get_position();
        const auto& velocity = object.get_velocity();
        m_output << i << "," << object.get_mass() << ",";
        m_output << position.get_X() << "," << position.get_Y() << "," << position.get_Z() << ",";
        m_output << velocity.get_X() << "," << velocity.get_Y() << "," << velocity.get_Z();
        m_output << "\n";
    }
}