#include "Constants.h"
#include "FastMultipoleMethod.h"

#include <algorithm>
#include <cmath>

using namespace nbs;

Vector3D FMMAlgorithm::convert_to_polar(const Vector3D& vec) {
    double rho = vec.get_modulus();
    if (rho < Constants::MINIMUM_RESOLUTION) {
        return Vector3D(0, 0, 0);
    }
    double alpha = std::acos(vec.get_Z() / rho);
    if (std::abs(vec.get_X()) < Constants::MINIMUM_RESOLUTION) {
        if (vec.get_Y() > 0) {
            return Vector3D(rho, alpha, Constants::PI / 2);
        } else {
            return Vector3D(rho, alpha, -Constants::PI / 2);
        }
    }
    double beta = std::atan2(vec.get_Y(), vec.get_X());
    return Vector3D(rho, alpha, beta);
}

Vector3D FMMAlgorithm::convert_to_cartesian(double A, double B, double C, const Vector3D& pos) {
    double alpha = pos.get_Y();
    double beta = pos.get_Z();
    double x = A * std::sin(alpha) * std::cos(beta) + B * std::cos(alpha) * std::cos(beta) -
               C * std::sin(beta);
    double y = A * std::sin(alpha) * std::sin(beta) + B * std::cos(alpha) * std::sin(beta) +
               C * std::cos(beta);
    double z = A * std::cos(alpha) - B * std::sin(alpha);
    return Vector3D(x, y, z);
}

bool FMMAlgorithm::is_adjacent(std::unique_ptr<OctreeNode<FMMData>>& node1,
                               std::unique_ptr<OctreeNode<FMMData>>& node2) {
    auto vec = node1->get_centre() - node2->get_centre();
    auto dx = std::abs(vec.get_X());
    auto dy = std::abs(vec.get_Y());
    auto dz = std::abs(vec.get_Z());
    auto x_range =
        node1->get_range().get_X() + node2->get_range().get_X() + Constants::MINIMUM_RESOLUTION;
    auto y_range =
        node1->get_range().get_Y() + node2->get_range().get_Y() + Constants::MINIMUM_RESOLUTION;
    auto z_range =
        node1->get_range().get_Z() + node2->get_range().get_Z() + Constants::MINIMUM_RESOLUTION;
    return dx < x_range && dy < y_range && dz < z_range;
}

void FMMAlgorithm::generate_L1_L3(std::unique_ptr<OctreeNode<FMMData>>& pos,    //
                                  std::unique_ptr<OctreeNode<FMMData>>& target, //
                                  std::vector<FMMNodeRef>& L1,                  //
                                  std::vector<FMMNodeRef>& L3) {
    if (pos->is_leaf()) {
        if (is_adjacent(pos, target)) {
            L1.emplace_back(pos);
            pos->get_data().get_list_L1().emplace_back(target);
        } else {
            L3.emplace_back(pos);
            pos->get_data().get_list_L4().emplace_back(target);
        }
    } else {
        if (!is_adjacent(pos, target)) {
            L3.emplace_back(pos);
            pos->get_data().get_list_L4().emplace_back(target);
        } else {
            for (const auto& child : OctreeRegion::all_children) {
                if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                    generate_L1_L3(ptr, target, L1, L3);
                }
            }
        }
    }
}

void FMMAlgorithm::generate_interaction_lists(std::unique_ptr<OctreeNode<FMMData>>& pos,
                                              OctreeNode<FMMData>* fa,
                                              const std::vector<FMMNodeRef>& fa_colleague) {
    auto& L2 = pos->get_data().get_list_L2();
    std::vector<FMMNodeRef> colleagues;
    if (fa != nullptr) {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = fa->get_child(child); ptr != nullptr && ptr != pos) {
                colleagues.emplace_back(ptr);
            }
        }
    }
    for (const auto& colleague : fa_colleague) {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = colleague.get()->get_child(child); ptr != nullptr) {
                if (is_adjacent(pos, ptr)) {
                    colleagues.emplace_back(ptr);
                } else {
                    L2.emplace_back(ptr);
                }
            }
        }
    }
    if (pos->is_leaf()) {
        auto& L1 = pos->get_data().get_list_L1();
        auto& L3 = pos->get_data().get_list_L3();
        for (const auto& colleague : colleagues) {
            if (colleague.get()->is_leaf()) {
                L1.emplace_back(colleague);
            } else {
                for (const auto& child : OctreeRegion::all_children) {
                    if (auto& ptr = colleague.get()->get_child(child); ptr != nullptr) {
                        generate_L1_L3(ptr, pos, L1, L3);
                    }
                }
            }
        }
    } else {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                generate_interaction_lists(ptr, pos.get(), colleagues);
            }
        }
    }
}

void FMMAlgorithm::calculate_P(double x, int p, std::vector<std::vector<double>>& Pnm) {
    Pnm.clear();
    Pnm.reserve(p);
    for (int n = 0; n < p; n++) {
        Pnm.emplace_back(std::vector<double>(n + 1, 0));
    }
    int dfac = 1;
    for (int m = 0; m < p; m++) {
        Pnm[m][m] = (m % 2 == 0) ? 1 : -1;
        Pnm[m][m] *= dfac;
        Pnm[m][m] *= std::pow((1 - x * x), m / 2.0);
        dfac *= (2 * m + 1);
    }
    for (int m = 0; m < p - 1; m++) {
        Pnm[m + 1][m] = (2 * m + 1) * x * Pnm[m][m];
    }
    for (int m = 0; m < p - 2; m++) {
        for (int n = m + 2; n < p; n++) {
            Pnm[n][m] = (2 * n - 1) * x * Pnm[n - 1][m] - (n + m - 1) * Pnm[n - 2][m];
            Pnm[n][m] /= (n - m);
        }
    }
}

void FMMAlgorithm::calculate_F(const Vector3D& vec, int p, std::vector<std::vector<Complex>>& F) {
    auto polar = convert_to_polar(vec);
    double rho = polar.get_X();
    double alpha = polar.get_Y();
    double beta = polar.get_Z();

    std::vector<std::vector<double>> P;
    calculate_P(std::cos(alpha), p, P);

    F.clear();
    F.reserve(p);
    for (int n = 0; n < p; n++) {
        F.emplace_back(std::vector<Complex>(n + 1, Complex(0, 0)));
    }
    for (int n = 0; n < p; n++) {
        for (int m = 0; m <= n; m++) {
            F[n][m] = Complex((n - m) % 2 == 0 ? 1 : -1, 0);
            F[n][m] *= std::pow(rho, n);
            F[n][m] /= Constants::FACTORIALS[n + m];
            F[n][m] *= std::exp(Complex(0, m * beta));
            F[n][m] *= P[n][m];
        }
    }
}

void FMMAlgorithm::calculate_G(const Vector3D& vec, int p, std::vector<std::vector<Complex>>& G) {
    auto polar = convert_to_polar(vec);
    double rho = polar.get_X();
    double alpha = polar.get_Y();
    double beta = polar.get_Z();

    std::vector<std::vector<double>> P;
    calculate_P(std::cos(alpha), p, P);

    G.clear();
    G.reserve(p);
    for (int n = 0; n < p; n++) {
        G.emplace_back(std::vector<Complex>(n + 1, Complex(0, 0)));
    }
    for (int n = 0; n < p; n++) {
        for (int m = 0; m <= n; m++) {
            G[n][m] = Complex((n - m) % 2 == 0 ? 1 : -1, 0);
            G[n][m] *= Constants::FACTORIALS[n - m];
            G[n][m] /= std::pow(rho, n + 1);
            G[n][m] *= std::exp(Complex(0, m * beta));
            G[n][m] *= P[n][m];
        }
    }
}

void FMMAlgorithm::init_expansion_vectors(std::unique_ptr<OctreeNode<FMMData>>& pos) {
    auto& M = pos->get_data().get_M();
    auto& L = pos->get_data().get_L();
    M.clear();
    L.clear();
    M.reserve(m_p);
    L.reserve(m_p);
    for (int n = 0; n < m_p; n++) {
        M.emplace_back(std::vector<Complex>(n + 1, Complex(0, 0)));
        L.emplace_back(std::vector<Complex>(n + 1, Complex(0, 0)));
    }
    for (const auto& child : OctreeRegion::all_children) {
        if (auto& ptr = pos->get_child(child); ptr != nullptr) {
            init_expansion_vectors(ptr);
        }
    }
}

void FMMAlgorithm::P2M_M2M(std::unique_ptr<OctreeNode<FMMData>>& pos) {
    auto& M = pos->get_data().get_M();
    if (pos->is_leaf()) {
        auto& objects = pos->get_objects();
        auto position = objects[0].get().get_position() - pos->get_centre();
        double total_mass = 0;
        for (const auto& object : objects) {
            total_mass += object.get().get_mass();
        }

        std::vector<std::vector<Complex>> F;
        calculate_F(-position, m_p, F);
        for (int n = 0; n < m_p; n++) {
            for (int m = 0; m <= n; m++) {
                M[n][m] = Complex(total_mass, 0);
                M[n][m] *= std::conj(F[n][m]);
            }
        }
    } else {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                P2M_M2M(ptr);
                auto vec_M2M = pos->get_centre() - ptr->get_centre();
                auto& child_M = ptr->get_data().get_M();
                std::vector<std::vector<Complex>> F;
                calculate_F(vec_M2M, m_p, F);
                for (int n1 = 0; n1 < m_p; n1++) {
                    for (int m1 = 0; m1 <= n1; m1++) {
                        for (int n = 0; n <= n1; n++) {
                            int m_min = std::max(-n, m1 - n1 + n);
                            int m_max = std::min(n, m1 + n1 - n);
                            for (int m = m_min; m <= m_max; m++) {
                                auto M_v =
                                    m >= 0 ? child_M[n][m]
                                           : std::conj(child_M[n][-m]) * (m % 2 == 0 ? 1.0 : -1.0);
                                auto F_v =
                                    m1 >= m ? std::conj(F[n1 - n][m1 - m])
                                            : F[n1 - n][m - m1] * ((m1 - m) % 2 == 0 ? 1.0 : -1.0);
                                M[n1][m1] += M_v * F_v;
                            }
                        }
                    }
                }
            }
        }
    }
}

Vector3D FMMAlgorithm::direct_calculation(PhysicalObject& target, PhysicalObject& source) {
    auto relative_position_vector = source.get_position() - target.get_position();
    double modulus = relative_position_vector.get_modulus();
    double partial_acceleration = Constants::G * source.get_mass() / (modulus * modulus * modulus);
    return relative_position_vector * partial_acceleration;
}

void FMMAlgorithm::P2P_P2L(std::unique_ptr<OctreeNode<FMMData>>& pos) {
    auto& L4 = pos->get_data().get_list_L4();
    if (auto& objects = pos->get_objects(); objects.size() <= m_p * m_p) {
        for (auto& object : objects) {
            Vector3D acceleration(0, 0, 0);
            for (auto& node : L4) {
                for (auto& source : node.get()->get_objects()) {
                    acceleration += direct_calculation(object.get(), source.get());
                }
            }
            const auto& acc_orig = object.get().get_acceleration();
            object.get().set_acceleration(acc_orig + acceleration);
        }
    } else {
        auto& L = pos->get_data().get_L();
        for (auto& node : L4) {
            auto& objects = node.get()->get_objects();
            auto position = objects[0].get().get_position() - pos->get_centre();
            double total_mass = 0;
            for (const auto& object : objects) {
                total_mass += object.get().get_mass();
            }

            std::vector<std::vector<Complex>> G;
            calculate_G(-position, m_p, G);
            for (int n = 0; n < m_p; n++) {
                for (int m = 0; m <= n; m++) {
                    L[n][m] += (Complex(total_mass, 0) * G[n][m]);
                }
            }
        }
    }
    for (const auto& child : OctreeRegion::all_children) {
        if (auto& ptr = pos->get_child(child); ptr != nullptr) {
            P2P_P2L(ptr);
        }
    }
}

void FMMAlgorithm::M2L(std::unique_ptr<OctreeNode<FMMData>>& pos) {
    auto& L = pos->get_data().get_L();
    auto& interaction_list = pos->get_data().get_list_L2();
    for (auto& node : interaction_list) {
        auto vec_M2L = pos->get_centre() - node.get()->get_centre();
        auto& M = node.get()->get_data().get_M();
        std::vector<std::vector<Complex>> G;
        calculate_G(vec_M2L, 2 * m_p - 1, G);
        for (int n1 = 0; n1 < m_p; n1++) {
            for (int m1 = 0; m1 <= n1; m1++) {
                for (int n = 0; n < m_p; n++) {
                    for (int m = -n; m <= n; m++) {
                        auto G_v = m1 + m >= 0 ? G[n1 + n][m1 + m]
                                               : std::conj(G[n1 + n][-m1 - m]) *
                                                     ((m1 + m) % 2 == 0 ? 1.0 : -1.0);
                        auto M_v =
                            m >= 0 ? M[n][m] : std::conj(M[n][-m]) * (m % 2 == 0 ? 1.0 : -1.0);
                        L[n1][m1] += G_v * M_v;
                    }
                }
            }
        }
    }
    for (const auto& child : OctreeRegion::all_children) {
        if (auto& ptr = pos->get_child(child); ptr != nullptr) {
            M2L(ptr);
        }
    }
}

void FMMAlgorithm::L2L(std::unique_ptr<OctreeNode<FMMData>>& pos, OctreeNode<FMMData>* fa) {
    if (fa != nullptr) {
        auto& L = pos->get_data().get_L();
        auto& parent_L = fa->get_data().get_L();
        auto vec_L2L = pos->get_centre() - fa->get_centre();
        std::vector<std::vector<Complex>> F;
        calculate_F(vec_L2L, m_p, F);
        for (int n1 = 0; n1 < m_p; n1++) {
            for (int m1 = 0; m1 <= n1; m1++) {
                for (int n = n1; n < m_p; n++) {
                    for (int m = m1 + n1 - n; m <= m1 + n - n1; m++) {
                        auto F_v = m >= m1 ? std::conj(F[n - n1][m - m1])
                                           : F[n - n1][m1 - m] * ((m - m1) % 2 == 0 ? 1.0 : -1.0);
                        auto L_v = m >= 0 ? parent_L[n][m]
                                          : std::conj(parent_L[n][-m]) * (m % 2 == 0 ? 1.0 : -1.0);
                        L[n1][m1] += F_v * L_v;
                    }
                }
            }
        }
    }
    for (const auto& child : OctreeRegion::all_children) {
        if (auto& ptr = pos->get_child(child); ptr != nullptr) {
            L2L(ptr, pos.get());
        }
    }
}

void FMMAlgorithm::P2P_L2P_M2P(std::unique_ptr<OctreeNode<FMMData>>& pos) {
    if (!pos->is_leaf()) {
        for (const auto& child : OctreeRegion::all_children) {
            if (auto& ptr = pos->get_child(child); ptr != nullptr) {
                P2P_L2P_M2P(ptr);
            }
        }
        return;
    }

    // P2P
    Vector3D acceleration(0, 0, 0);
    auto& objects = pos->get_objects();
    auto& list_L1 = pos->get_data().get_list_L1();
    for (auto& node : list_L1) {
        for (auto& source : node.get()->get_objects()) {
            acceleration += direct_calculation(objects[0].get(), source.get());
        }
    }

    // L2P
    auto L2P_vec = objects[0].get().get_position() - pos->get_centre();
    auto polar = convert_to_polar(L2P_vec);
    double rho = polar.get_X();
    double alpha = polar.get_Y();
    double beta = polar.get_Z();
    auto& L = pos->get_data().get_L();
    if (rho < Constants::MINIMUM_RESOLUTION) {
        acceleration -= Constants::G * Vector3D(-L[1][1].real(), -L[1][1].imag(), L[1][0].real());
    } else if (alpha < Constants::MINIMUM_RESOLUTION) {
        Vector3D acc(0, 0, 0);
        double fac = -1;
        for (int n = 1; n < m_p; n++) {
            acc += Vector3D(fac * L[n][1].real(), fac * L[n][1].imag(), fac * L[n][0].real());
            fac *= -rho / n;
        }
        acceleration -= Constants::G * acc;
    } else if (Constants::PI - alpha < Constants::MINIMUM_RESOLUTION) {
        Vector3D acc(0, 0, 0);
        double fac = -1;
        for (int n = 1; n < m_p; n++) {
            acc += Vector3D(fac * L[n][1].real(), fac * L[n][1].imag(), -fac * L[n][0].real());
            fac *= -rho / n;
        }
        acceleration -= Constants::G * acc;
    } else {
        std::vector<std::vector<Complex>> F;
        calculate_F(L2P_vec, m_p, F);
        double A_tot = 0;
        double B_tot = 0;
        double C_tot = 0;
        for (int n = 1; n < m_p; n++) {
            A_tot += (L[n][0] * F[n][0]).real() * n / rho;
            B_tot -=
                (std::exp(Complex(0, beta)) * std::conj(F[n][1]) * L[n][0]).real() * (n + 1) / rho;
            for (int m = 1; m <= n; m++) {
                auto grad_F = F[n][m] * (n / rho);
                A_tot += 2 * L[n][m].real() * grad_F.real() + 2 * L[n][m].imag() * grad_F.imag();
                grad_F = F[n][m] * (m / std::tan(alpha));
                if (m < n) {
                    grad_F -= (n + m + 1) * 1.0 * std::exp(Complex(0, -beta)) * F[n][m + 1];
                }
                grad_F /= rho;
                B_tot += 2 * L[n][m].real() * grad_F.real() + 2 * L[n][m].imag() * grad_F.imag();
                grad_F = Complex(0, m) * F[n][m] / rho / std::sin(alpha);
                C_tot += 2 * L[n][m].real() * grad_F.real() + 2 * L[n][m].imag() * grad_F.imag();
            }
        }
        acceleration -= Constants::G * convert_to_cartesian(A_tot, B_tot, C_tot, polar);
    }

    // M2P
    auto& L3 = pos->get_data().get_list_L3();
    for (auto& node : L3) {
        if (auto& sources = node.get()->get_objects(); sources.size() <= m_p * m_p) {
            for (auto& source : sources) {
                acceleration += direct_calculation(objects[0].get(), source.get());
            }
        } else {
            auto M2P_vec = objects[0].get().get_position() - pos->get_centre();
            polar = convert_to_polar(M2P_vec);
            rho = polar.get_X();
            alpha = polar.get_Y();
            beta = polar.get_Z();
            auto& M = node.get()->get_data().get_M();
            if (alpha < Constants::MINIMUM_RESOLUTION) {
                double fac = -1.0 / rho / rho;
                Vector3D acc(0, 0, fac * M[0][0].real());
                for (int n = 1; n < m_p; n++) {
                    acc += Vector3D(-fac * M[n][1].real(), fac * M[n][1].imag(),
                                    -fac * M[n][0].real());
                    fac *= (-(n + 1) / rho);
                }
                acceleration += Constants::G * acc;
            } else if (Constants::PI - alpha < Constants::MINIMUM_RESOLUTION) {
                double fac = -1.0 / rho / rho;
                Vector3D acc(0, 0, -fac * M[0][0].real());
                for (int n = 1; n < m_p; n++) {
                    acc +=
                        Vector3D(-fac * M[n][1].real(), fac * M[n][1].imag(), fac * M[n][0].real());
                    fac *= (-(n + 1) / rho);
                }
                acceleration += Constants::G * acc;
            } else {
                std::vector<std::vector<Complex>> G;
                calculate_G(M2P_vec, m_p, G);
                double A_tot = -(M[0][0] * G[0][0]).real() / rho;
                double B_tot = 0;
                double C_tot = 0;
                for (int n = 1; n < m_p; n++) {
                    A_tot -= (M[n][0] * G[n][0]).real() * (n + 1) / rho;
                    B_tot -= (std::exp(Complex(0, -beta)) * G[n][1] * M[n][0]).real() * n / rho;
                    for (int m = 1; m <= n; m++) {
                        auto grad_G = G[n][m] * (-(n + 1) / rho);
                        A_tot +=
                            2 * M[n][m].real() * grad_G.real() - 2 * M[n][m].imag() * grad_G.imag();
                        grad_G = G[n][m] * (m / std::tan(alpha));
                        if (m < n) {
                            grad_G -= (n - m) * 1.0 * std::exp(Complex(0, -beta)) * G[n][m + 1];
                        }
                        grad_G /= rho;
                        B_tot +=
                            2 * M[n][m].real() * grad_G.real() - 2 * M[n][m].imag() * grad_G.imag();
                        grad_G = Complex(0, m) * G[n][m] / rho / std::sin(alpha);
                        C_tot +=
                            2 * M[n][m].real() * grad_G.real() - 2 * M[n][m].imag() * grad_G.imag();
                    }
                }
                acceleration += Constants::G * convert_to_cartesian(A_tot, B_tot, C_tot, polar);
            }
        }
    }

    // Update acceleration
    for (const auto& object : objects) {
        const auto& acc_orig = object.get().get_acceleration();
        object.get().set_acceleration(acc_orig + acceleration);
    }
}

void FastMultipoleMethod::calculate_accelerations(std::vector<PhysicalObject>& objects) {
    m_logger.log("Accelerations calculation started.");
    auto root = std::make_unique<OctreeNode<FMMData>>( //
        Vector3D(m_x_range, m_y_range, m_z_range), Vector3D(0, 0, 0));
    std::vector<Object_ref> object_refs;
    for (auto& object : objects) {
        object.set_acceleration(Vector3D(0, 0, 0));
        object_refs.emplace_back(object);
    }
    root->set_objects(std::move(object_refs));
    Octree<FMMData> octree(std::move(root));
    m_logger.log("Octree building started.");
    octree.build_octree();
    m_logger.log("Octree building finished.");

    FMMAlgorithm algorithm(octree, m_p);
    m_logger.log("FMM algorithm started.");
    m_logger.log("Initialization of expansion vectors started.");
    algorithm.init_expansion_vectors();
    m_logger.log("Initialization of expansion vectors finished.");
    m_logger.log("Generation of interaction lists started.");
    algorithm.generate_interaction_lists();
    m_logger.log("Generation of interaction lists finished.");
    m_logger.log("P2M and M2M started.");
    algorithm.P2M_M2M();
    m_logger.log("P2M and M2M finished.");
    m_logger.log("P2P and P2L started.");
    algorithm.P2P_P2L();
    m_logger.log("P2P and P2L finished.");
    m_logger.log("M2L started.");
    algorithm.M2L();
    m_logger.log("M2L finished.");
    m_logger.log("L2L started.");
    algorithm.L2L();
    m_logger.log("L2L finished.");
    m_logger.log("P2P, L2P and M2P started.");
    algorithm.P2P_L2P_M2P();
    m_logger.log("P2P, L2P and M2P finished.");
    m_logger.log("FMM algorithm finished.");
}