#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/reference/ReferencePlatform.h"
#include "ReferenceMSFKernelFactory.h"
#include "ReferenceMSFKernels.h"

using namespace MSFPlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if(platform.getName() == "Reference") {
//        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceMSFKernelFactory* factory = new ReferenceMSFKernelFactory();
            platform.registerKernelFactory(CalcShrinkForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerMSFReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* ReferenceMSFKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& referencePlatformData = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
//    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());

    if (name == CalcShrinkForceKernel::Name())
        return new ReferenceCalcShrinkForceKernel(name, platform);

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str());
}