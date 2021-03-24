#ifdef WIN32
#define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/ShrinkForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "ReferenceMSFKernelFactory.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/ShrinkKernels.h"
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace MSFPlugin;
using namespace std;

const double TOL = 1e-5;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceMSFKernelFactory* factory = new ReferenceMSFKernelFactory();
            platform.registerKernelFactory(CalcShrinkForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerMSFReferenceKernelFactories() {
    registerKernelFactories();
}

void testSimpleExpression() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    ShrinkForce* forceField = new ShrinkForce("-0.1*r^3");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    system.addForce(forceField);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
//    vector<Vec3> positions(2);
//    positions[0] = Vec3(0, 0, 0);
//    positions[1] = Vec3(2, 0, 0);
//    context.setPositions(positions);
//    State state = context.getState(State::Forces | State::Energy);
//    const vector<Vec3>& forces = state.getForces();
//    double force = 0.1*3*(2*2);
//    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
//    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
//    ASSERT_EQUAL_TOL(-0.1*(2*2*2), state.getPotentialEnergy(), TOL);
}

int main() {
    try {
        registerMSFReferenceKernelFactories();
        testSimpleExpression();
//        testParameter();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}