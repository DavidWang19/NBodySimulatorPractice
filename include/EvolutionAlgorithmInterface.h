#pragma once

#include "OutputUtils.h"
#include "PhysicalObject.h"

namespace nbs
{
    class EvolutionAlgorithmInterface
    {
    public:
        EvolutionAlgorithmInterface() = delete;
        EvolutionAlgorithmInterface(EvolutionAlgorithmInterface&&) = default;
        ~EvolutionAlgorithmInterface() = default;
        EvolutionAlgorithmInterface(const std::string& logger) : m_logger(logger) {}

        void initialize(std::vector<PhysicalObject>& objects, double time_step);
        void evolve(std::vector<PhysicalObject>& objects, double time_step);

    protected:
        virtual void other_initializations() {}
        virtual void calculate_accelerations(std::vector<PhysicalObject>& objects) = 0;

        Logger m_logger;
    };
}