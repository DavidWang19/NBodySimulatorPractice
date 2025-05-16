#include "EvolutionAlgorithmInterface.h"

#include <format>

using namespace nbs;

void EvolutionAlgorithmInterface::initialize(std::vector<PhysicalObject>& objects,
                                             double time_step) {
    m_logger.log("Initialization started.");
    calculate_accelerations(objects);
    for (auto& object : objects) {
        object.set_velocity(object.get_velocity() + object.get_acceleration() * time_step * 0.5);
    }
    other_initializations();
    m_logger.log("Initialization finished.");
}

void EvolutionAlgorithmInterface::evolve(std::vector<PhysicalObject>& objects, double time_step) {
    m_logger.log(std::format("Evolution for t = {}s started.", time_step));
    for (auto& object : objects) {
        object.set_position(object.get_position() + object.get_velocity() * time_step);
    }
    calculate_accelerations(objects);
    for (auto& object : objects) {
        object.set_velocity(object.get_velocity() + object.get_acceleration() * time_step);
    }
    m_logger.log(std::format("Evolution for t = {}s finished.", time_step));
}