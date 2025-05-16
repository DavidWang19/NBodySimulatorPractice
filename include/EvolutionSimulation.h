#pragma once

#include "EvolutionAlgorithmInterface.h"
#include "PhysicalObject.h"

namespace nbs
{
    class EvolutionSimulation
    {
    public:
        EvolutionSimulation() = default;
        EvolutionSimulation(EvolutionSimulation&&) = default;
        ~EvolutionSimulation() = default;
        EvolutionSimulation(std::vector<PhysicalObject>&& objects,
                            std::unique_ptr<EvolutionAlgorithmInterface>&& evolution_algorithm_ptr,
                            double time_step)
            : m_objects(std::move(objects)),
              m_evolution_algorithm_ptr(std::move(evolution_algorithm_ptr)),
              m_time_step(time_step) {
            m_evolution_algorithm_ptr->initialize(m_objects, m_time_step);
        };

        double get_time() const { return m_time; }
        const std::vector<PhysicalObject>& get_objects() const { return m_objects; }

        void evolve() {
            m_evolution_algorithm_ptr->evolve(m_objects, m_time_step);
            m_time += m_time_step;
        }

    protected:
        std::vector<PhysicalObject> m_objects;
        std::unique_ptr<EvolutionAlgorithmInterface> m_evolution_algorithm_ptr;
        double m_time_step;
        double m_time = 0;
    };
}