#pragma once

#include "EvolutionAlgorithmInterface.h"

namespace nbs
{
    class NaiveMethod : public EvolutionAlgorithmInterface
    {
    public:
        NaiveMethod() = delete;
        NaiveMethod(NaiveMethod&&) = default;
        ~NaiveMethod() = default;
        NaiveMethod(const std::string& logger) : EvolutionAlgorithmInterface(logger) {}

    private:
        void calculate_accelerations(std::vector<PhysicalObject>& objects) override;
    };
}
