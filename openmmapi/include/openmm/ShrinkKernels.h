#ifndef OPENMM_SHRINKKERNELS_H_
#define OPENMM_SHRINKKERNELS_H_

#include "openmm/ShrinkForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include <iosfwd>
#include <set>
#include <string>
#include <vector>


namespace MSFPlugin {
/**
 * This kernel is invoked by ShrinkForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcShrinkForceKernel : public OpenMM::KernelImpl {
    public:
        enum NonbondedMethod {
            NoCutoff = 0,
            CutoffNonPeriodic = 1,
            CutoffPeriodic = 2,
        };

        static std::string Name() {
            return "CalcShrinkForce";
        }

        CalcShrinkForceKernel(std::string name, const OpenMM::Platform &platform) : KernelImpl(name, platform) {
        }

        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the ShrinkForce this kernel will be used for
         */
        virtual void initialize(const OpenMM::System &system, const ShrinkForce &force) = 0;

        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        virtual double execute(OpenMM::ContextImpl &context, bool includeForces, bool includeEnergy) = 0;

        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the ShrinkForce to copy the parameters from
         */
        virtual void copyParametersToContext(OpenMM::ContextImpl &context, const ShrinkForce &force) = 0;
    };
} // namespace MSFPlugin

#endif /*OPENMM_SHRINKKERNELS_H_*/