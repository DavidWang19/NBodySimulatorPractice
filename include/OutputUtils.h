#pragma once

#include "Octree.h"

#include <chrono>
#include <fstream>
#include <string_view>

namespace nbs
{
    class OutputUtils
    {
    public:
        OutputUtils() = delete;
        OutputUtils(OutputUtils&&) = delete;
        ~OutputUtils();
        OutputUtils(const std::string& file_name);

    protected:
        std::ofstream m_output;
    };

    class Logger : public OutputUtils
    {
    public:
        Logger() = delete;
        Logger(Logger&&) = delete;
        ~Logger() = default;
        Logger(const std::string& file_name) : OutputUtils(file_name) {
            m_last_time = std::chrono::system_clock::now();
            log("Logger started");
        };

        void log(std::string_view message);

    private:
        std::chrono::system_clock::time_point m_last_time;
    };

    template <typename T>
    class OctreeOutput
    {
    public:
        OctreeOutput() = delete;
        ~OctreeOutput() { m_output.close(); }
        OctreeOutput(Octree<T>& target, int mode, const std::string& file_name) : m_target(target) {
            if (mode == 0) {
                m_output.open(file_name);
            } else if (mode == 1) {
                m_output.open(file_name, std::ios::app);
            }
            m_output.setf(std::ios::fixed, std::ios::floatfield);
            m_output.precision(4);
        };

        void output_octree() {
            output_node(m_target.get_root());
            m_output << "********************************\n";
        }

    private:
        std::ofstream m_output;
        Octree<T>& m_target;

        void output_node(const std::unique_ptr<OctreeNode<T>>& node) {
            auto x = node->get_centre().get_X() - node->get_range().get_X();
            auto y = node->get_centre().get_Y() - node->get_range().get_Y();
            auto w = node->get_range().get_X() * 2;
            auto h = node->get_range().get_Y() * 2;
            m_output << x << "," << y << "," << w << "," << h << "\n";
            for (const auto& child : OctreeRegion::all_children) {
                if (auto& ptr = node->get_child(child); ptr != nullptr) {
                    output_node(ptr);
                }
            }
        }
    };
} // namespace nbs