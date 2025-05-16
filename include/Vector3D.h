#pragma once

namespace nbs
{
    class Vector3D
    {
    public:
        Vector3D() = default;
        ~Vector3D() = default;
        Vector3D(const Vector3D&) = default;
        Vector3D(Vector3D&&) = default;
        Vector3D(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

        double get_X() const { return m_x; }
        double get_Y() const { return m_y; }
        double get_Z() const { return m_z; }
        void set_X(double x) { m_x = x; }
        void set_Y(double y) { m_y = y; }
        void set_Z(double z) { m_z = z; }

        void operator=(const Vector3D& rhs) {
            m_x = rhs.m_x;
            m_y = rhs.m_y;
            m_z = rhs.m_z;
        }

        void operator=(Vector3D&& rhs) {
            m_x = rhs.m_x;
            m_y = rhs.m_y;
            m_z = rhs.m_z;
        }

        Vector3D& operator+=(const Vector3D& rhs);
        Vector3D& operator-=(const Vector3D& rhs);
        Vector3D& operator*=(double scalar);
        friend Vector3D operator+(Vector3D lhs, const Vector3D& rhs);
        friend Vector3D operator-(Vector3D lhs, const Vector3D& rhs);
        friend Vector3D operator-(Vector3D vec);
        friend Vector3D operator*(Vector3D lhs, double scalar);
        friend Vector3D operator*(double scalar, Vector3D rhs);

        double dot(const Vector3D& rhs);
        Vector3D cross(const Vector3D& rhs);

        double get_modulus() const;

    protected:
        double m_x;
        double m_y;
        double m_z;
    };
} // namespace nbs
