#include "Constants.h"
#include "NaiveMethod.h"
#include "Vector3D.h"

using namespace nbs;

void NaiveMethod::calculate_accelerations(std::vector<PhysicalObject>& objects) {
    for (auto& object : objects) {
        object.set_acceleration(Vector3D(0, 0, 0));
    }
    for (size_t i = 0; i < objects.size(); i++) {
        for (size_t j = i + 1; j < objects.size(); j++) {
            Vector3D relative_position_vector =
                objects[j].get_position() - objects[i].get_position();
            double modulus = relative_position_vector.get_modulus();
            double factor = Constants::G / (modulus * modulus * modulus);
            objects[i].set_acceleration(objects[i].get_acceleration() +
                                        relative_position_vector * factor * objects[j].get_mass());
            objects[j].set_acceleration(objects[j].get_acceleration() -
                                        relative_position_vector * factor * objects[i].get_mass());
        }
    }
}