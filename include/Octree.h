#pragma once

#include "Constants.h"
#include "PhysicalObject.h"

#include <array>
#include <functional>
#include <initializer_list>
#include <memory>
#include <utility>
#include <vector>

namespace nbs
{
    using Object_ref = std::reference_wrapper<PhysicalObject>;

    class OctreeRegion
    {
    public:
        enum Value
        {
            PX_PY_PZ, // 000 = 0
            PX_PY_NZ, // 001 = 1
            PX_NY_PZ, // 010 = 2
            PX_NY_NZ, // 011 = 3
            NX_PY_PZ, // 100 = 4
            NX_PY_NZ, // 101 = 5
            NX_NY_PZ, // 110 = 6
            NX_NY_NZ  // 111 = 7
        };

        static Value get_region(const Vector3D& rel_vec);
        static std::pair<Vector3D, Vector3D> get_new_configuration //
            (Value region, const Vector3D& parent_range, const Vector3D& parent_centre);

        static constexpr std::initializer_list<Value> all_children {
            PX_PY_PZ, PX_PY_NZ, PX_NY_PZ, PX_NY_NZ, //
            NX_PY_PZ, NX_PY_NZ, NX_NY_PZ, NX_NY_NZ  //
        };
    };

    template <typename T>
    class OctreeNode
    {
    public:
        OctreeNode() = delete;
        ~OctreeNode() = default;
        OctreeNode(Vector3D&& range, Vector3D&& centre)
            : m_range(std::move(range)), m_centre(std::move(centre)) {};

        std::unique_ptr<OctreeNode<T>>& get_child(OctreeRegion::Value region) {
            return m_children[static_cast<int>(region)];
        }

        std::vector<Object_ref>& get_objects() { return m_objects; }
        void set_objects(std::vector<Object_ref>&& objects) { m_objects = std::move(objects); }
        Vector3D& get_range() { return m_range; }
        Vector3D& get_centre() { return m_centre; }
        void set_is_leaf(bool is_leaf) { m_is_leaf = is_leaf; }
        bool is_leaf() { return m_is_leaf; }
        T& get_data() { return m_data; }

    protected:
        std::array<std::unique_ptr<OctreeNode<T>>, 8> m_children;
        std::vector<Object_ref> m_objects;
        Vector3D m_range;
        Vector3D m_centre;
        bool m_is_leaf = false;
        T m_data;
    };

    template <typename T>
    class Octree
    {
    public:
        Octree() = delete;
        ~Octree() = default;
        Octree(std::unique_ptr<OctreeNode<T>>&& root) : m_root(std::move(root)) {};

        std::unique_ptr<OctreeNode<T>>& get_root() { return m_root; }

        void build_octree() { build_octree(m_root); }

    protected:
        std::unique_ptr<OctreeNode<T>> m_root;

        void build_octree(std::unique_ptr<OctreeNode<T>>& pos) {
            auto& objects = pos->get_objects();
            if (objects.size() == 1) {
                pos->set_is_leaf(true);
                return;
            }
            if (pos->get_range().get_X() < Constants::MINIMUM_RESOLUTION) {
                pos->set_is_leaf(true);
                return;
            }
            auto& centre = pos->get_centre();
            for (auto& object : objects) {
                auto region = OctreeRegion::get_region(object.get().get_position() - centre);
                auto& child_ptr = pos->get_child(region);
                if (child_ptr == nullptr) {
                    auto [new_range, new_centre] = OctreeRegion::get_new_configuration(
                        region, pos->get_range(), pos->get_centre());
                    std::unique_ptr<OctreeNode<T>> new_child = std::make_unique<OctreeNode<T>>(
                        std::move(new_range), std::move(new_centre));
                    child_ptr = std::move(new_child);
                }
                child_ptr->get_objects().emplace_back(object);
            }
            for (const auto& child : OctreeRegion::all_children) {
                if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                    build_octree(ptr);
                }
            }
        }
    };

} // namespace nbs