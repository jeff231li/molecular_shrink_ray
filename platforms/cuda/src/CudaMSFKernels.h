#ifndef OPENMM_CUDAMSFKERNELS_H_
#define OPENMM_CUDAMSFKERNELS_H_

#include "CudaPlatform.h"
#include "CudaContext.h"
#include "CudaKernels.h"
#include "openmm/common/CommonKernels.h"


namespace OpenMM {
/**
 * This kernel is invoked by ShrinkForce to calculate the forces acting on the system.
 */
    class CudaCalcShrinkForceKernel : public CalcShrinkForceKernel {
    public:
        CudaParallelCalcShrinkForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data,
                                     const System& system);
        CommonCalcShrinkForceKernel& getKernel(int index) {
            return dynamic_cast<CommonCalcShrinkForceKernel&>(kernels[index].getImpl());
        }
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
        class Task;
        CudaPlatform::PlatformData& data;
        std::vector<Kernel> kernels;
    };

} // namespace OpenMM

#endif /*OPENMM_CUDAMSFKERNELS_H_*/