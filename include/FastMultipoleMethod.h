#pragma once

#include "EvolutionAlgorithmInterface.h"
#include "Octree.h"

#include <complex>
#include <vector>

namespace nbs
{
    using Complex = std::complex<double>;

    class FMMData
    {
    public:
        FMMData() = default;
        ~FMMData() = default;

        std::vector<std::vector<Complex>>& get_M() { return m_M; }
        std::vector<std::vector<Complex>>& get_L() { return m_L; }
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>>& get_list_L1() {
            return m_list_L1;
        }
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>>& get_list_L2() {
            return m_list_L2;
        }
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>>& get_list_L3() {
            return m_list_L3;
        }
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>>& get_list_L4() {
            return m_list_L4;
        }

    protected:
        std::vector<std::vector<Complex>> m_M;
        std::vector<std::vector<Complex>> m_L;
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>> m_list_L1;
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>> m_list_L2;
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>> m_list_L3;
        std::vector<std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>> m_list_L4;
    };

    using FMMNodeRef = std::reference_wrapper<std::unique_ptr<OctreeNode<FMMData>>>;

    class FMMAlgorithm
    {
    public:
        FMMAlgorithm() = delete;
        ~FMMAlgorithm() = default;
        FMMAlgorithm(Octree<FMMData>& target, int p) : m_target(target), m_p(p) {}

        void init_expansion_vectors() { init_expansion_vectors(m_target.get_root()); }
        void generate_interaction_lists() {
            generate_interaction_lists(m_target.get_root(), nullptr, std::vector<FMMNodeRef>());
        }

        void P2M_M2M() { P2M_M2M(m_target.get_root()); }
        void P2P_P2L() { P2P_P2L(m_target.get_root()); }
        void M2L() { M2L(m_target.get_root()); }
        void L2L() { L2L(m_target.get_root(), nullptr); }
        void P2P_L2P_M2P() { P2P_L2P_M2P(m_target.get_root()); }

    private:
        Octree<FMMData>& m_target;
        int m_p;

        Vector3D convert_to_polar(const Vector3D& vec);
        Vector3D convert_to_cartesian(double A, double B, double C, const Vector3D& pos);
        bool is_adjacent(std::unique_ptr<OctreeNode<FMMData>>& node1,
                         std::unique_ptr<OctreeNode<FMMData>>& node2);
        void calculate_P(double x, int p, std::vector<std::vector<double>>& Pnm);
        void calculate_F(const Vector3D& vec, int p, std::vector<std::vector<Complex>>& F);
        void calculate_G(const Vector3D& vec, int p, std::vector<std::vector<Complex>>& G);
        void init_expansion_vectors(std::unique_ptr<OctreeNode<FMMData>>& pos);
        void generate_L1_L3(std::unique_ptr<OctreeNode<FMMData>>& pos,    //
                            std::unique_ptr<OctreeNode<FMMData>>& target, //
                            std::vector<FMMNodeRef>& L1,                  //
                            std::vector<FMMNodeRef>& L3);
        void generate_interaction_lists(std::unique_ptr<OctreeNode<FMMData>>& pos,
                                        OctreeNode<FMMData>* fa,
                                        const std::vector<FMMNodeRef>& fa_colleague);
        Vector3D direct_calculation(PhysicalObject& target, PhysicalObject& source);
        void P2M_M2M(std::unique_ptr<OctreeNode<FMMData>>& pos);
        void P2P_P2L(std::unique_ptr<OctreeNode<FMMData>>& pos);
        void M2L(std::unique_ptr<OctreeNode<FMMData>>& pos);
        void L2L(std::unique_ptr<OctreeNode<FMMData>>& pos, OctreeNode<FMMData>* fa);
        void P2P_L2P_M2P(std::unique_ptr<OctreeNode<FMMData>>& pos);
    };

    class FastMultipoleMethod : public EvolutionAlgorithmInterface
    {
    public:
        FastMultipoleMethod() = delete;
        FastMultipoleMethod(FastMultipoleMethod&&) = default;
        ~FastMultipoleMethod() = default;
        FastMultipoleMethod(double x_range, double y_range, double z_range, int p,
                            const std::string& logger)
            : m_x_range(x_range), m_y_range(y_range), m_z_range(z_range), m_p(p),
              EvolutionAlgorithmInterface(logger) {}

    private:
        double m_x_range;
        double m_y_range;
        double m_z_range;
        int m_p;

        void calculate_accelerations(std::vector<PhysicalObject>& objects) override;
    };
} // namespace nbs