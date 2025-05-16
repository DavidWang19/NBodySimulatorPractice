#include "Vector3D.h"

#include <cmath>

using namespace nbs;

Vector3D& Vector3D::operator+=(const Vector3D& rhs) {
    m_x += rhs.m_x;
    m_y += rhs.m_y;
    m_z += rhs.m_z;
    return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D& rhs) {
    m_x -= rhs.m_x;
    m_y -= rhs.m_y;
    m_z -= rhs.m_z;
    return *this;
}

Vector3D& Vector3D::operator*=(double scalar) {
    m_x *= scalar;
    m_y *= scalar;
    m_z *= scalar;
    return *this;
}

Vector3D nbs::operator+(Vector3D lhs, const Vector3D& rhs) {
    lhs += rhs;
    return lhs;
}

Vector3D nbs::operator-(Vector3D lhs, const Vector3D& rhs) {
    lhs -= rhs;
    return lhs;
}

Vector3D nbs::operator-(Vector3D vec) {
    vec *= -1;
    return vec;
}

Vector3D nbs::operator*(Vector3D lhs, double scalar) {
    lhs *= scalar;
    return lhs;
}

Vector3D nbs::operator*(double scalar, Vector3D rhs) {
    rhs *= scalar;
    return rhs;
}

double Vector3D::dot(const Vector3D& rhs) {
    return m_x * rhs.m_x + m_y * rhs.m_y + m_z * rhs.m_z;
}

Vector3D Vector3D::cross(const Vector3D& rhs) {
    return Vector3D(m_y * rhs.m_z - m_z * rhs.m_y, //
                    m_z * rhs.m_x - m_x * rhs.m_z, //
                    m_x * rhs.m_y - m_y * rhs.m_x);
}

double Vector3D::get_modulus() const {
    return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
}