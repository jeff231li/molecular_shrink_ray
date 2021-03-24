#ifndef OPENMM_CPUMSFKERNELS_H_
#define OPENMM_CPUMSFKERNELS_H_

#include "CpuShrinkForce.h"
#include "openmm/CpuPlatform.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include <array>
#include <tuple>

namespace OpenMM {
/**
 * This kernel is invoked by ShrinkForce to calculate the forces acting on the system.
 */
    class CpuCalcShrinkForceKernel : public CalcShrinkForceKernel {
    public:
        CpuCalcShrinkForceKernel(std::string name, const Platform &platform, CpuPlatform::PlatformData &data);

        ~CpuCalcShrinkForceKernel();

        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the ShrinkForce this kernel will be used for
         */
        void initialize(const System &system, const ShrinkForce &force);

        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        double execute(ContextImpl &context, bool includeForces, bool includeEnergy);

        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the ShrinkForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl &context, const ShrinkForce &force);

    private:
        CpuPlatform::PlatformData &data;
        int numParticles;
        std::vector <std::vector<double>> particleParamArray;
        double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
        bool useSwitchingFunction, hasInitializedLongRangeCorrection;
        ShrinkForce *forceCopy;
        std::map<std::string, double> globalParamValues;
        std::vector <std::set<int>> exclusions;
        std::vector <std::string> parameterNames, globalParameterNames, energyParamDerivNames;
        std::vector <std::pair<std::set < int>, std::set<int>> >
        interactionGroups;
        std::vector<double> longRangeCoefficientDerivs;
        NonbondedMethod nonbondedMethod;
        CpuShrinkForce *nonbonded;
    };
} // namespace OpenMM

#endif /*OPENMM_CPUMSFKERNELS_H_*/