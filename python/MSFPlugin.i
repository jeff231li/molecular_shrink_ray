%module MSFPlugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */
%include "std_vector.i"
%include "std_set.i"
namespace std {
    %template(vectord) vector<double>;
    %template(vectori) vector<int>;
    %template(seti) set<int>;
}

%{
#include "OpenMMMSF.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
 */
%pythonprepend MSFPlugin::ShrinkForce::setCutoffDistance(double distance) %{
    if isinstance(distance, unit.Quantity):
        distance = distance.value_in_unit(unit.nanometer)
%}
%pythonappend MSFPlugin::ShrinkForce::getCutoffDistance() const %{
    val = unit.Quantity(val, unit.nanometer)
%}
%pythonprepend MSFPlugin::ShrinkForce::setSwitchingDistance(double distance) %{
    if isinstance(distance, unit.Quantity):
        distance = distance.value_in_unit(unit.nanometer)
%}
%pythonappend MSFPlugin::ShrinkForce::getSwitchingDistance() const %{
    val = unit.Quantity(val, unit.nanometer)
%}
//%pythonprepend MSFPlugin::ShrinkForce::addParticle(const std::vector<double>& parameters=std::vector<double>()) %{
//    for parm in args:
//        if isinstance(parm, unit):
//%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

namespace MSFPlugin{

class ShrinkForce : public OpenMM::Force {
    public:
        enum NonbondedMethod {
            NoCutoff = 0,
            CutoffNonPeriodic = 1,
            CutoffPeriodic = 2,
        };
        ShrinkForce(const std::string& energy);

        int getNumParticles() const;
        int getNumExclusions() const;
        int getNumPerParticleParameters() const;
        int getNumGlobalParameters() const;
        int getNumTabulatedFunctions() const;
        int getNumFunctions() const;
        int getNumInteractionGroups() const;
        int getNumEnergyParameterDerivatives() const;
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
        %apply std::vector<double>& OUTPUT {std::vector<double>& parameters};
        void getParticleParameters(int index, std::vector<double>& parameters) const;
        %clear std::vector<double>& parameters;
        void setParticleParameters(int index, const std::vector<double>& parameters);

        int addExclusion(int particle1, int particle2);
        %apply int& OUTPUT {int& particle1};
        %apply int& OUTPUT {int& particle2};
        void getExclusionParticles(int index, int& particle1, int& particle2) const;
        %clear int& particle1;
        %clear int& particle2;
        void setExclusionParticles(int index, int particle1, int particle2);
        void createExclusionsFromBonds(const std::vector<std::pair<int, int> >& bonds, int bondCutoff);

        int addTabulatedFunction(const std::string& name, TabulatedFunction* function);
        const TabulatedFunction& getTabulatedFunction(int index) const;
        TabulatedFunction& getTabulatedFunction(int index);
        const std::string& getTabulatedFunctionName(int index) const;

        int addFunction(const std::string& name, const std::vector<double>& values, double min, double max);
        %apply std::string& OUTPUT {std::string& name};
        %apply std::vector<double>& OUTPUT {std::vector<double>& values};
        %apply double& OUTPUT {double& min};
        %apply double& OUTPUT {double& max};
        void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const;
        %clear std::string& name;
        %clear std::vector<double>& values;
        %clear double& min;
        %clear double& max;
        void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max);

        int addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2);
        %apply std::set<int>& OUTPUT {std::set<int>& set1};
        %apply std::set<int>& OUTPUT {std::set<int>& set2};
        void getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const;
        %clear std::set<int>& set1;
        %clear std::set<int>& set2;
        void setInteractionGroupParameters(int index, const std::set<int>& set1, std::set<int>& set2);

        void updateParametersInContext(Context& context);
        bool usesPeriodicBoundaryConditions() const;

        /*
         * Add methods for casting a Force to a TorchForce.
         */
        %extend {
            static MSFPlugin::ShrinkForce& cast(OpenMM::Force& force) {
                return dynamic_cast<MSFPlugin::ShrinkForce&>(force);
            }

            static bool isinstance(OpenMM::Force& force) {
                return (dynamic_cast<MSFPlugin::ShrinkForce*>(&force) != NULL);
            }
        }
    };
}