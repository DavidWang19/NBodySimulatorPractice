#pragma once

#include "Vector3D.h"

#include <memory>

namespace nbs
{
    class PhysicalObject
    {
    public:
        PhysicalObject() = default;
        PhysicalObject(PhysicalObject&&) = default;
        ~PhysicalObject() = default;
        PhysicalObject(double mass, Vector3D&& position, Vector3D&& velocity)
            : m_mass(mass), m_position(std::move(position)), m_velocity(std::move(velocity)),
              m_acceleration(0, 0, 0) {};
        PhysicalObject(double mass, double x, double y, double z, double vx, double vy, double vz)
            : m_mass(mass), m_position(x, y, z), m_velocity(vx, vy, vz), m_acceleration(0, 0, 0) {};

        bool operator==(const PhysicalObject& other) const { return &other == this; }

        double get_mass() const { return m_mass; }
        const Vector3D& get_position() const { return m_position; }
        const Vector3D& get_velocity() const { return m_velocity; }
        const Vector3D& get_acceleration() const { return m_acceleration; }
        void set_position(const Vector3D& position) { m_position = position; }
        void set_velocity(const Vector3D& velocity) { m_velocity = velocity; }
        void set_acceleration(const Vector3D& acceleration) { m_acceleration = acceleration; }

    protected:
        double m_mass;
        Vector3D m_position;
        Vector3D m_velocity;
        Vector3D m_acceleration;
    };
}
