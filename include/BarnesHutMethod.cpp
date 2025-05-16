#include "BarnesHutMethod.h"

#include "Constants.h"
#include "OutputUtils.h"

using namespace nbs;

void BHAlgorithm::calculate_mass_centre() {
    calculate_mass_centre(m_target.get_root());
}

void BHAlgorithm::calculate_mass_centre(std::unique_ptr<OctreeNode<BHData>>& pos) {
    auto& objects = pos->get_objects();
    double total_mass = 0;
    Vector3D mass_centre(0, 0, 0);
    for (auto& object : objects) {
        total_mass += object.get().get_mass();
        mass_centre += object.get().get_position() * object.get().get_mass();
    }
    mass_centre *= (1.0 / total_mass);
    pos->get_data().set_total_mass(total_mass);
    pos->get_data().set_mass_centre(mass_centre);
    for (const auto& child : OctreeRegion::all_children) {
        if (auto& ptr = pos->get_child(child); ptr != nullptr) {
            calculate_mass_centre(ptr);
        }
    }
}

void BarnesHutMethod::other_initializations() {
    this->m_output_mode = 1;
}

void BarnesHutMethod::calculate_acceleration( //
    std::unique_ptr<OctreeNode<BHData>>& pos, //
    PhysicalObject& object,                   //
    Vector3D& result                          //
) {
    auto& centre = pos->get_centre();
    auto& mass_centre = pos->get_data().get_mass_centre();
    auto& range = pos->get_range();
    auto relative_position_vector = mass_centre - object.get_position();
    double modulus = relative_position_vector.get_modulus();
    if (modulus < Constants::MINIMUM_RESOLUTION) return;
    if (pos->is_leaf() || range.get_modulus() / modulus < m_theta) {
        double partial_acceleration =
            Constants::G * pos->get_data().get_total_mass() / (modulus * modulus * modulus);
        result += relative_position_vector * partial_acceleration;
    } else {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                calculate_acceleration(ptr, object, result);
            }
        }
    }
}

void BarnesHutMethod::calculate_accelerations(std::vector<PhysicalObject>& objects) {
    m_logger.log("Accelerations calculation started.");
    auto root = std::make_unique<OctreeNode<BHData>>( //
        Vector3D(m_x_range, m_y_range, m_z_range), Vector3D(0, 0, 0));
    std::vector<Object_ref> object_refs;
    for (auto& object : objects) {
        object_refs.emplace_back(object);
    }
    root->set_objects(std::move(object_refs));
    Octree<BHData> octree(std::move(root));
    m_logger.log("Octree building started.");
    octree.build_octree();
    m_logger.log("Octree building finished.");
    BHAlgorithm algorithm(octree);
    m_logger.log("Mass centre calculation started.");
    algorithm.calculate_mass_centre();
    m_logger.log("Mass centre calculation finished.");
    m_logger.log("Individual accelerations calculation started.");
    for (auto& object : objects) {
        Vector3D acceleration(0, 0, 0);
        calculate_acceleration(octree.get_root(), object, acceleration);
        object.set_acceleration(acceleration);
    }
    m_logger.log("Individual accelerations calculation finished.");

    if (m_output_file_name == "") return;
    OctreeOutput output(octree, m_output_mode, m_output_file_name);
    output.output_octree();
}