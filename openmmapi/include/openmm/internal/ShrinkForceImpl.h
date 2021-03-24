#ifndef OPENMM_SHRINKFORCEIMPL_H_
#define OPENMM_SHRINKFORCEIMPL_H_

#include "openmm/internal/ForceImpl.h"
#include "openmm/ShrinkForce.h"
#include "openmm/ShrinkKernels.h"
#include "openmm/Kernel.h"
#include "lepton/CompiledExpression.h"
#include <utility>
#include <map>
#include <string>


namespace MSFPlugin {
/**
 * This is the internal implementation of ShrinkForce.
 */
class OPENMM_EXPORT ShrinkForceImpl : public OpenMM::ForceImpl {
public:
    ShrinkForceImpl(const ShrinkForce& owner);
    ~ShrinkForceImpl();
    void initialize(OpenMM::ContextImpl& context);
    const ShrinkForce& getOwner() const {
        return owner;
    }
    void updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) {}
    double calcForcesAndEnergy(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(OpenMM::ContextImpl& context);

    static void calcLongRangeCorrection(const ShrinkForce& force, const OpenMM::Context& context, double& coefficient,
                                        std::vector<double>& derivatives);
private:
    static double integrateInteraction(Lepton::CompiledExpression& expression, const std::vector<double>& params1,
                                       const std::vector<double>& params2, const ShrinkForce& force, const OpenMM::Context&
                                       context, const std::vector<std::string>& paramNames);
    const ShrinkForce& owner;
    OpenMM::Kernel kernel;
};
} // namespace MSFPlugin

#endif /*OPENMM_SHRINKFORCEIMPL_H_*/