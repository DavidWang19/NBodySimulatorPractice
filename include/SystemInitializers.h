#pragma once

#include "PhysicalObject.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace nbs
{
    class ReadFromFile
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, const std::string& filename) {
            std::ifstream file(filename);
            std::string line;
            while (std::getline(file, line)) {
                if (!std::isdigit(line[0])) break;
                auto comma1 = line.find(',');
                auto comma2 = line.find(',', comma1 + 1);
                double mass = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                comma2 = line.find(',', comma1 + 1);
                double x = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                comma2 = line.find(',', comma1 + 1);
                double y = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                comma2 = line.find(',', comma1 + 1);
                double z = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                comma2 = line.find(',', comma1 + 1);
                double vx = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                comma2 = line.find(',', comma1 + 1);
                double vy = std::stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
                comma1 = comma2;
                double vz = std::stod(line.substr(comma1 + 1));
                objects.emplace_back(mass, x, y, z, vx, vy, vz);
            }
            file.close();
        }
    };

    class RandomDistribution
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, int num_objects,
                               double mass_mean, double mass_stddev, double velocity_mean,
                               double velocity_stddev, double x_range, double y_range,
                               double z_range) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> mass_distribution(mass_mean, mass_stddev);
            std::normal_distribution<double> velocity_distribution(velocity_mean, velocity_stddev);
            std::normal_distribution<double> standard_distribution(0, 1);
            std::uniform_real_distribution<double> x_distribution(-x_range, x_range);
            std::uniform_real_distribution<double> y_distribution(-y_range, y_range);
            std::uniform_real_distribution<double> z_distribution(-z_range, z_range);
            for (int i = 0; i < num_objects; ++i) {
                auto velocity = velocity_distribution(gen);
                auto vx = standard_distribution(gen);
                auto vy = standard_distribution(gen);
                auto vz = standard_distribution(gen);
                auto mod = std::sqrt(vx * vx + vy * vy + vz * vz);
                vx /= mod;
                vy /= mod;
                vz /= mod;
                vx *= velocity;
                vy *= velocity;
                vz *= velocity;
                double z = z_range == 0 ? 0 : z_distribution(gen);
                double mass = std::abs(mass_distribution(gen));
                objects.emplace_back(mass, x_distribution(gen), y_distribution(gen), z, vx, vy, vz);
            }
        }
    };

    class UniformStatic2D
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, int num_objects, double mass,
                               double x_range, double y_range) {
            int objects_per_side = std::floor(std::sqrt(num_objects));
            double x_step = 2 * x_range / (objects_per_side + 1);
            double y_step = 2 * y_range / (objects_per_side + 1);
            for (int i = 0; i * i < num_objects; ++i) {
                for (int j = 0; j * j < num_objects; ++j) {
                    objects.emplace_back(mass, (i + 0.5) * x_step - x_range,
                                         (j + 0.5) * y_step - y_range, 0, 0, 0, 0);
                }
            }
        }
    };

    class UniformStatic3D
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, int num_objects, double mass,
                               double x_range, double y_range, double z_range) {
            int objects_per_side = std::floor(std::cbrt(num_objects));
            double x_step = 2 * x_range / (objects_per_side + 1);
            double y_step = 2 * y_range / (objects_per_side + 1);
            double z_step = 2 * z_range / (objects_per_side + 1);
            for (int i = 0; i * i * i < num_objects; ++i) {
                for (int j = 0; j * j * j < num_objects; ++j) {
                    for (int k = 0; k * k * k < num_objects; ++k) {
                        objects.emplace_back(mass, (i + 0.5) * x_step - x_range,
                                             (j + 0.5) * y_step - y_range,
                                             (k + 0.5) * z_step - z_range, 0, 0, 0);
                    }
                }
            }
        }
    };

    class Inhomogeneous2D
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, int num_objects,
                               int num_clusters, double range) {
            int orig_num_objects = num_objects;
            std::vector<std::pair<nbs::Vector3D, int>> cluster_centers;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> pos_distribution(-range, range);
            for (int i = 0; i < num_clusters - 1; ++i) {
                std::uniform_int_distribution<int> num_objects_distribution(1, num_objects);
                int num_objects_in_cluster = num_objects_distribution(gen);
                num_objects -= num_objects_in_cluster;
                cluster_centers.emplace_back(
                    nbs::Vector3D(pos_distribution(gen), pos_distribution(gen), 0),
                    num_objects_in_cluster);
            }
            cluster_centers.emplace_back(
                nbs::Vector3D(pos_distribution(gen), pos_distribution(gen), 0), num_objects);
            for (const auto& cluster : cluster_centers) {
                std::normal_distribution<double> mass_distribution(1, 0.1);
                std::normal_distribution<double> velocity_distribution(0, 0.1);
                std::normal_distribution<double> standard_distribution(0, 1);
                double coverage = std::sqrt(cluster.second * 1.0 / orig_num_objects);
                coverage *= range * 0.1;
                std::normal_distribution<double> dev_from_center_distribution(-coverage, coverage);
                for (int i = 0; i < cluster.second; ++i) {
                    auto velocity = velocity_distribution(gen);
                    auto vx = standard_distribution(gen);
                    auto vy = standard_distribution(gen);
                    auto mod = std::sqrt(vx * vx + vy * vy);
                    vx /= mod;
                    vy /= mod;
                    vx *= velocity;
                    vy *= velocity;
                    double mass = std::abs(mass_distribution(gen));
                    objects.emplace_back(
                        mass, cluster.first.get_X() + dev_from_center_distribution(gen),
                        cluster.first.get_Y() + dev_from_center_distribution(gen), 0, vx, vy, 0);
                }
            }
        }
    };

    class Inhomogeneous3D
    {
    public:
        static void initialize(std::vector<PhysicalObject>& objects, int num_objects,
                               int num_clusters, double range) {
            int orig_num_objects = num_objects;
            std::vector<std::pair<nbs::Vector3D, int>> cluster_centers;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> pos_distribution(-range, range);
            for (int i = 0; i < num_clusters - 1; ++i) {
                std::uniform_int_distribution<int> num_objects_distribution(1, num_objects);
                int num_objects_in_cluster = num_objects_distribution(gen);
                num_objects -= num_objects_in_cluster;
                cluster_centers.emplace_back(nbs::Vector3D(pos_distribution(gen),
                                                           pos_distribution(gen),
                                                           pos_distribution(gen)),
                                             num_objects_in_cluster);
            }
            cluster_centers.emplace_back(
                nbs::Vector3D(pos_distribution(gen), pos_distribution(gen), pos_distribution(gen)),
                num_objects);
            for (const auto& cluster : cluster_centers) {
                std::normal_distribution<double> mass_distribution(1, 0.1);
                std::normal_distribution<double> velocity_distribution(0, 0.1);
                std::normal_distribution<double> standard_distribution(0, 1);
                double coverage = std::cbrt(cluster.second * 1.0 / orig_num_objects);
                coverage *= range * 0.1;
                std::normal_distribution<double> dev_from_center_distribution(-coverage, coverage);
                for (int i = 0; i < cluster.second; ++i) {
                    auto velocity = velocity_distribution(gen);
                    auto vx = standard_distribution(gen);
                    auto vy = standard_distribution(gen);
                    auto vz = standard_distribution(gen);
                    auto mod = std::sqrt(vx * vx + vy * vy + vz * vz);
                    vx /= mod;
                    vy /= mod;
                    vz /= mod;
                    vx *= velocity;
                    vy *= velocity;
                    vz *= velocity;
                    double mass = std::abs(mass_distribution(gen));
                    objects.emplace_back(
                        mass, cluster.first.get_X() + dev_from_center_distribution(gen),
                        cluster.first.get_Y() + dev_from_center_distribution(gen),
                        cluster.first.get_Z() + dev_from_center_distribution(gen), vx, vy, vz);
                }
            }
        }
    };
} // namespace nbs