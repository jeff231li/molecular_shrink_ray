#ifndef OPENMM_SHRINKFORCE_H_
#define OPENMM_SHRINKFORCE_H_

#include "openmm/TabulatedFunction.h"
#include "openmm/Force.h"
#include "openmm/Vec3.h"
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "internal/windowsExportMSF.h"

using namespace OpenMM;

namespace MSFPlugin {
/**
 * This class is a fork of CustomNonbondedForce that includes centroid-based calculations
 * like in CustomCentroidBondForce. This can be used to "Shrink" a molecule.
 */
class OPENMM_EXPORT ShrinkForce : public OpenMM::Force {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2,
    };

    // Constructor and destructor
    explicit ShrinkForce(const std::string& energy);
    ShrinkForce(const ShrinkForce& rhs);
    ~ShrinkForce();

    int getNumParticles() const {
        return particles.size();
    }
    int getNumExclusions() const {
        return exclusions.size();
    }
    int getNumPerParticleParameters() const {
        return parameters.size();
    }
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    int getNumTabulatedFunctions() const {
        return functions.size();
    }
    int getNumFunctions() const {
        return functions.size();
    }
    int getNumInteractionGroups() const {
        return interactionGroups.size();
    }
    int getNumEnergyParameterDerivatives() const {
        return energyParameterDerivatives.size();
    }
    const std::string& getEnergyFunction() const;
    void setEnergyFunction(const std::string& energy);

    NonbondedMethod getNonbondedMethod() const;
    void setNonbondedMethod(NonbondedMethod method);

    double getCutoffDistance() const;
    void setCutoffDistance(double distance);

    bool getUseSwitchingFunction() const;
    void setUseSwitchingFunction(bool use);

    double getSwitchingDistance() const;
    void setSwitchingDistance(double distance);

    bool getUseLongRangeCorrection() const;
    void setUseLongRangeCorrection(bool use);

    int addPerParticleParameter(const std::string& name);
    const std::string& getPerParticleParameterName(int index) const;
    void setPerParticleParameterName(int index, const std::string& name);

    int addGlobalParameter(const std::string& name, double defaultValue);
    const std::string& getGlobalParameterName(int index) const;
    void setGlobalParameterName(int index, const std::string& name);
    double getGlobalParameterDefaultValue(int index) const;
    void setGlobalParameterDefaultValue(int index, double defaultValue);

    void addEnergyParameterDerivative(const std::string& name);
    const std::string& getEnergyParameterDerivativeName(int index) const;

    int addParticle(const std::vector<double>& parameters=std::vector<double>());
    void getParticleParameters(int index, std::vector<double>& parameters) const;
    void setParticleParameters(int index, const std::vector<double>& parameters);

    int addExclusion(int particle1, int particle2);
    void getExclusionParticles(int index, int& particle1, int& particle2) const;
    void setExclusionParticles(int index, int particle1, int particle2);

    void createExclusionsFromBonds(const std::vector<std::pair<int, int> >& bonds, int bondCutoff);

    int addTabulatedFunction(const std::string& name, TabulatedFunction* function);
    const TabulatedFunction& getTabulatedFunction(int index) const;
    TabulatedFunction& getTabulatedFunction(int index);
    const std::string& getTabulatedFunctionName(int index) const;

    int addFunction(const std::string& name, const std::vector<double>& values, double min, double max);
    void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const;
    void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max);

    int addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2);
    void getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const;
    void setInteractionGroupParameters(int index, const std::set<int>& set1, std::set<int>& set2);

    void updateParametersInContext(OpenMM::Context& context);

    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == ShrinkForce::CutoffPeriodic;
    }
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class PerParticleParameterInfo;
    class GlobalParameterInfo;
    class ExclusionInfo;
    class FunctionInfo;
    class InteractionGroupInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, switchingDistance;
    bool useSwitchingFunction, useLongRangeCorrection;
    std::string energyExpression;
    std::vector<PerParticleParameterInfo> parameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<ParticleInfo> particles;
    std::vector<ExclusionInfo> exclusions;
    std::vector<FunctionInfo> functions;
    std::vector<InteractionGroupInfo> interactionGroups;
    std::vector<int> energyParameterDerivatives;
};

class ShrinkForce::ParticleInfo {
public:
    std::vector<double> parameters;
    ParticleInfo() {}
    ParticleInfo(const std::vector<double>& parameters) : parameters(parameters) {}
};

class ShrinkForce::PerParticleParameterInfo {
public:
    std::string name;
    PerParticleParameterInfo() {}
    PerParticleParameterInfo(const std::string& name) : name(name) {}
};

class ShrinkForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {}
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {}
};

class ShrinkForce::ExclusionInfo {
public:
    int particle1, particle2;
    ExclusionInfo() {
        particle1 = particle2 = -1;
    }
    ExclusionInfo(int particle1, int particle2) :
        particle1(particle1), particle2(particle2) {}
};

class ShrinkForce::FunctionInfo {
public:
    std::string name;
    OpenMM::TabulatedFunction* function;
    FunctionInfo() {}
    FunctionInfo(const std::string& name, OpenMM::TabulatedFunction* function) : name(name), function(function) {}
};

class ShrinkForce::InteractionGroupInfo {
public:
    std::set<int> set1, set2;
    InteractionGroupInfo() {}
    InteractionGroupInfo(const std::set<int>& set1, const std::set<int>& set2) : set1(set1), set2(set2) {}
};

} // namespace MSFPlugin

#endif /*OPENMM_SHRINKFORCE_H_*/
