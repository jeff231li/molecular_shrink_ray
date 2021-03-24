#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/ShrinkForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ShrinkForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace MSFPlugin;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

ShrinkForce::ShrinkForce(const string& energy) : energyExpression(energy), nonbondedMethod(NoCutoff),
cutoffDistance(1.0), switchingDistance(-1.0), useSwitchingFunction(false), useLongRangeCorrection(false) {}

ShrinkForce::ShrinkForce(const ShrinkForce& rhs) {
    energyExpression = rhs.energyExpression;
    nonbondedMethod = rhs.nonbondedMethod;
    cutoffDistance = rhs.cutoffDistance;
    switchingDistance = rhs.switchingDistance;
    useSwitchingFunction = rhs.useSwitchingFunction;
    useLongRangeCorrection = rhs.useLongRangeCorrection;
    parameters = rhs.parameters;
    globalParameters = rhs.globalParameters;
    energyParameterDerivatives = rhs.energyParameterDerivatives;
    particles = rhs.particles;
    exclusions = rhs.exclusions;
    interactionGroups = rhs.interactionGroups;
    for (std::vector<FunctionInfo>::const_iterator it = rhs.functions.begin(); it != rhs.functions.end(); ++it)
        functions.push_back(FunctionInfo(it->name, it->function->Copy()));
}
ShrinkForce::~ShrinkForce() {
    for (auto function : functions)
        delete function.function;
}

const string& ShrinkForce::getEnergyFunction() const {
    return energyExpression;
}
void ShrinkForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

ShrinkForce::NonbondedMethod ShrinkForce::getNonbondedMethod() const {
    return nonbondedMethod;
}
void ShrinkForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("ShrinkForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double ShrinkForce::getCutoffDistance() const {
    return cutoffDistance;
}
void ShrinkForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

bool ShrinkForce::getUseSwitchingFunction() const {
    return useSwitchingFunction;
}
void ShrinkForce::setUseSwitchingFunction(bool use) {
    useSwitchingFunction = use;
}
double ShrinkForce::getSwitchingDistance() const {
    return switchingDistance;
}
void ShrinkForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}

bool ShrinkForce::getUseLongRangeCorrection() const {
    return useLongRangeCorrection;
}
void ShrinkForce::setUseLongRangeCorrection(bool use) {
    useLongRangeCorrection = use;
}

int ShrinkForce::addPerParticleParameter(const std::string& name) {
    parameters.push_back(PerParticleParameterInfo(name));
    return parameters.size()-1;
}
const std::string& ShrinkForce::getPerParticleParameterName(int index) const {
    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}
void ShrinkForce::setPerParticleParameterName(int index, const std::string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int ShrinkForce::addGlobalParameter(const std::string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}
const std::string& ShrinkForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}
void ShrinkForce::setGlobalParameterName(int index, const std::string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}
double ShrinkForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}
void ShrinkForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void ShrinkForce::addEnergyParameterDerivative(const std::string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(std::string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}
const std::string& ShrinkForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int ShrinkForce::addParticle(const std::vector<double>& parameters) {
    particles.push_back(ParticleInfo(parameters));
    return particles.size()-1;
}
void ShrinkForce::getParticleParameters(int index, std::vector<double>& parameters) const {
    ASSERT_VALID_INDEX(index, particles);
    parameters = particles[index].parameters;
}
void ShrinkForce::setParticleParameters(int index, const std::vector<double>& parameters) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].parameters = parameters;
}

int ShrinkForce::addExclusion(int particle1, int particle2) {
    exclusions.push_back(ExclusionInfo(particle1, particle2));
    return exclusions.size()-1;
}
void ShrinkForce::getExclusionParticles(int index, int& particle1, int& particle2) const {
    ASSERT_VALID_INDEX(index, exclusions);
    particle1 = exclusions[index].particle1;
    particle2 = exclusions[index].particle2;
}
void ShrinkForce::setExclusionParticles(int index, int particle1, int particle2) {
    ASSERT_VALID_INDEX(index, exclusions);
    exclusions[index].particle1 = particle1;
    exclusions[index].particle2 = particle2;
}
void ShrinkForce::createExclusionsFromBonds(const std::vector<std::pair<int, int>>& bonds, int bondCutoff) {
    if (bondCutoff < 1)
        return;
    for (auto& bond : bonds)
        if (bond.first < 0 || bond.second < 0 || bond.first >= particles.size() || bond.second >= particles.size())
            throw OpenMMException("createExclusionsFromBonds: Illegal particle index in list of bonds");
    std::vector<std::set<int> > exclusions(particles.size());
    std::vector<std::set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        int p1 = bond.first;
        int p2 = bond.second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bondCutoff-1; level++) {
        std::vector<std::set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) particles.size(); i++)
            for (int j : currentExclusions[i])
                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i)
                addExclusion(i, j);
}

int ShrinkForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}
const TabulatedFunction& ShrinkForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}
TabulatedFunction& ShrinkForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}
const string& ShrinkForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

int ShrinkForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max) {
    functions.push_back(FunctionInfo(name, new Continuous1DFunction(values, min, max)));
    return functions.size()-1;
}
void ShrinkForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min,
                                        double& max) const {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("ShrinkForce: function is not a Continuous1DFunction");
    name = functions[index].name;
    function->getFunctionParameters(values, min, max);
}
void ShrinkForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double
min, double max) {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("ShrinkForce: function is not a Continuous1DFunction");
    functions[index].name = name;
    function->setFunctionParameters(values, min, max);
}

int ShrinkForce::addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2) {
    for (std::set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
        ASSERT(*it >= 0);
    for (std::set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
        ASSERT(*it >= 0);
    interactionGroups.push_back(InteractionGroupInfo(set1, set2));
    return interactionGroups.size()-1;
}
void ShrinkForce::getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const {
    ASSERT_VALID_INDEX(index, interactionGroups);
    set1 = interactionGroups[index].set1;
    set2 = interactionGroups[index].set2;
}
void ShrinkForce::setInteractionGroupParameters(int index, const std::set<int>& set1, std::set<int>& set2) {
    ASSERT_VALID_INDEX(index, interactionGroups);
    for (std::set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
        ASSERT_VALID_INDEX(*it, particles);
    for (std::set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
        ASSERT_VALID_INDEX(*it, particles);
    interactionGroups[index].set1 = set1;
    interactionGroups[index].set2 = set2;
}

ForceImpl* ShrinkForce::createImpl() const {
    return new ShrinkForceImpl(*this);
}

void ShrinkForce::updateParametersInContext(OpenMM::Context& context) {
    dynamic_cast<ShrinkForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
