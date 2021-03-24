#ifndef __ReferenceShrinkxIxn_H__
#define __ReferenceShrinkxIxn_H__

#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/reference/ReferencePairIxn.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {

    class ReferenceShrinkIxn {

    private:

        bool cutoff;
        bool useSwitch;
        bool periodic;
        const OpenMM::NeighborList* neighborList;
        OpenMM::Vec3 periodicBoxVectors[3];
        double cutoffDistance, switchingDistance;
        Lepton::CompiledExpression energyExpression;
        Lepton::CompiledExpression forceExpression;
        std::vector<std::string> paramNames;
        std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
        CompiledExpressionSet expressionSet;
        std::vector<int> particleParamIndex;
        int rIndex;
        std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;

        /**---------------------------------------------------------------------------------------

           Calculate custom pair ixn between two atoms

           @param atom1            the index of the first atom
           @param atom2            the index of the second atom
           @param atomCoordinates  atom coordinates
           @param forces           force array (forces added)
           @param totalEnergy      total energy

           --------------------------------------------------------------------------------------- */

        void calculateOneIxn(int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces,
                             double* totalEnergy, double* energyParamDerivs);


    public:

        /**---------------------------------------------------------------------------------------

           Constructor

           --------------------------------------------------------------------------------------- */

        ReferenceShrinkIxn(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression&
        forceExpression, const std::vector<std::string>& parameterNames, const
                           std::vector<Lepton::CompiledExpression> energyParamDerivExpressions);

        /**---------------------------------------------------------------------------------------

           Destructor

           --------------------------------------------------------------------------------------- */

        ~ReferenceShrinkIxn();

        /**---------------------------------------------------------------------------------------

           Set the force to use a cutoff.

           @param distance            the cutoff distance
           @param neighbors           the neighbor list to use

           --------------------------------------------------------------------------------------- */

        void setUseCutoff(double distance, const OpenMM::NeighborList& neighbors);

        /**---------------------------------------------------------------------------------------

           Restrict the force to a list of interaction groups.

           @param distance            the cutoff distance
           @param neighbors           the neighbor list to use

           --------------------------------------------------------------------------------------- */

        void setInteractionGroups(const std::vector<std::pair<std::set<int>, std::set<int> > >& groups);

        /**---------------------------------------------------------------------------------------

           Set the force to use a switching function.

           @param distance            the switching distance

           --------------------------------------------------------------------------------------- */

        void setUseSwitchingFunction(double distance);

        /**---------------------------------------------------------------------------------------

           Set the force to use periodic boundary conditions.  This requires that a cutoff has
           already been set, and the smallest side of the periodic box is at least twice the cutoff
           distance.

           @param vectors    the vectors defining the periodic box

           --------------------------------------------------------------------------------------- */

        void setPeriodic(OpenMM::Vec3* vectors);

        /**---------------------------------------------------------------------------------------

           Calculate custom pair ixn

           @param numberOfAtoms    number of atoms
           @param atomCoordinates  atom coordinates
           @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
           @param exclusions       atom exclusion indices
                                   exclusions[atomIndex] contains the list of exclusions for that atom
           @param globalParameters the values of global parameters
           @param forces           force array (forces added)
           @param totalEnergy      total energy

           --------------------------------------------------------------------------------------- */

        void calculatePairIxn(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                              std::vector<std::vector<double> >& atomParameters, std::vector<std::set<int> >& exclusions,
                              const std::map<std::string, double>& globalParameters, std::vector<OpenMM::Vec3>& forces,
                              double* totalEnergy, double* energyParamDerivs);

// ---------------------------------------------------------------------------------------

    };

} // namespace OpenMM

#endif // __ReferenceShrinkxIxn_H__
