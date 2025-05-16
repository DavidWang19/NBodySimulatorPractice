#pragma once

#include "EvolutionSimulation.h"
#include "OutputUtils.h"

namespace nbs
{
    class SystemOutput : public OutputUtils
    {
    public:
        SystemOutput() = delete;
        SystemOutput(SystemOutput&&) = delete;
        ~SystemOutput() = default;
        SystemOutput(const std::string& file_name) : OutputUtils(file_name) {};

        void output_system_arrangement(EvolutionSimulation& simulation);
        void output_system_arrangement_raw(const std::vector<PhysicalObject>& objects);
    };
}