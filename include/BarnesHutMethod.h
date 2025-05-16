#pragma once

#include "EvolutionAlgorithmInterface.h"
#include "Octree.h"

#include <string>

namespace nbs
{
    class BHData
    {
    public:
        BHData() = default;
        ~BHData() = default;
        BHData(double total_mass, const Vector3D& mass_centre)
            : m_total_mass(total_mass), m_mass_centre(mass_centre) {};

        double get_total_mass() { return m_total_mass; }
        Vector3D& get_mass_centre() { return m_mass_centre; }
        void set_total_mass(double mass) { m_total_mass = mass; }
        void set_mass_centre(const Vector3D& centre) { m_mass_centre = centre; }

    protected:
        double m_total_mass;
        Vector3D m_mass_centre;
    };

    class BHAlgorithm
    {
    public:
        BHAlgorithm() = delete;
        ~BHAlgorithm() = default;
        BHAlgorithm(Octree<BHData>& target) : m_target(target) {};

        void calculate_mass_centre();

    private:
        Octree<BHData>& m_target;

        void calculate_mass_centre(std::unique_ptr<nbs::OctreeNode<nbs::BHData>>& pos);
    };

    class BarnesHutMethod : public EvolutionAlgorithmInterface
    {
    public:
        BarnesHutMethod() = delete;
        BarnesHutMethod(BarnesHutMethod&&) = default;
        ~BarnesHutMethod() = default;
        BarnesHutMethod(double x_range, double y_range, double z_range, double theta,
                        const std::string& logger)
            : m_x_range(x_range), m_y_range(y_range), m_z_range(z_range), m_theta(theta),
              EvolutionAlgorithmInterface(logger) {}

        void set_output_file_name(const std::string& output_file_name) {
            m_output_file_name = output_file_name;
        }

    private:
        double m_x_range;
        double m_y_range;
        double m_z_range;
        double m_theta;
        std::string m_output_file_name = "";
        int m_output_mode = 0; // 0: overwrite, 1: append

        void other_initializations() override;

        void calculate_acceleration(std::unique_ptr<OctreeNode<BHData>>& pos,
                                    PhysicalObject& object, Vector3D& result);
        void calculate_accelerations(std::vector<PhysicalObject>& objects) override;
    };
} // namespace nbs