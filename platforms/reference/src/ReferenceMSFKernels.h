#ifndef OPENMM_REFERENCEMSFKERNELS_H_
#define OPENMM_REFERENCEMSFKERNELS_H_

#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/ShrinkKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "openmm/reference/ReferencePlatform.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <array>
#include <utility>

using namespace OpenMM;
using namespace std;

namespace MSFPlugin {

/**
 * This kernel is invoked by ShrinkForce to calculate the forces acting on the system.
 */
class ReferenceCalcShrinkForceKernel : public CalcShrinkForceKernel {
public:
    ReferenceCalcShrinkForceKernel(std::string name, const Platform& platform) :
            CalcShrinkForceKernel(name, platform), forceCopy(NULL) {
    }
    ~ReferenceCalcShrinkForceKernel();

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the ShrinkForce this kernel will be used for
     */
    void initialize(const System& system, const ShrinkForce& force);

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);

    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ShrinkForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const ShrinkForce& force);

    private:
        int numParticles;
        std::vector<std::vector<double> > particleParamArray;
        double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
        bool useSwitchingFunction, hasInitializedLongRangeCorrection;
        ShrinkForce *forceCopy;
        std::map<std::string, double> globalParamValues;
        std::vector<std::set<int> > exclusions;
        Lepton::CompiledExpression energyExpression, forceExpression;
        std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
        std::vector<std::string> parameterNames, globalParameterNames, energyParamDerivNames;
        std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;
        std::vector<double> longRangeCoefficientDerivs;
        NonbondedMethod nonbondedMethod;
        NeighborList* neighborList;
    };

} // namespace MSFPlugin

#endif /*OPENMM_REFERENCEMSFKERNELS_H_*/