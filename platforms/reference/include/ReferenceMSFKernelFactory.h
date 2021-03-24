#ifndef OPENMM_REFERENCEMSFKERNELFACTORY_H_
#define OPENMM_REFERENCEMSFKERNELFACTORY_H_


#include "openmm/KernelFactory.h"

namespace OpenMM {

/**
 * This KernelFactory creates kernels for the reference implementation of the MSF plugin
 */

class ReferenceMSFKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEMSFKERNELFACTORY_H_*/