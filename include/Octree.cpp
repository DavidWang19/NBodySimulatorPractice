#include "Octree.h"

using namespace nbs;

OctreeRegion::Value OctreeRegion::get_region(const Vector3D& rel_vec) {
    int xf = rel_vec.get_X() < 0;
    int yf = rel_vec.get_Y() < 0;
    int zf = rel_vec.get_Z() < 0;
    return static_cast<OctreeRegion::Value>((xf << 2) + (yf << 1) + zf);
}

std::pair<Vector3D, Vector3D> OctreeRegion::get_new_configuration //
    (OctreeRegion::Value region, const Vector3D& parent_range, const Vector3D& parent_centre) {
    double new_x_r = parent_range.get_X() / 2;
    double new_y_r = parent_range.get_Y() / 2;
    double new_z_r = parent_range.get_Z() / 2;
    int region_int = static_cast<int>(region);
    int xf = 1 - (region_int >> 2) * 2;
    int yf = 1 - ((region_int >> 1) & 1) * 2;
    int zf = 1 - (region_int & 1) * 2;
    double new_x_c = parent_centre.get_X() + xf * new_x_r;
    double new_y_c = parent_centre.get_Y() + yf * new_y_r;
    double new_z_c = parent_centre.get_Z() + zf * new_z_r;
    return std::make_pair(Vector3D(new_x_r, new_y_r, new_z_r), Vector3D(new_x_c, new_y_c, new_z_c));
}