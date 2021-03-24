#include "CudaMSFKernels.h"
#include "CudaKernelSources.h"


using namespace OpenMM;
using namespace std;


class CudaCalcShrinkForceKernel::Task : public CudaContext::WorkTask {
public:
    Task(ContextImpl& context, CommonCalcShrinkForceKernel& kernel, bool includeForce,
         bool includeEnergy, double& energy) : context(context), kernel(kernel),
                                               includeForce(includeForce), includeEnergy(includeEnergy), energy(energy) {
    }
    void execute() {
        energy += kernel.execute(context, includeForce, includeEnergy);
    }
private:
    ContextImpl& context;
    CommonCalcShrinkForceKernel& kernel;
    bool includeForce, includeEnergy;
    double& energy;
};

CudaCalcShrinkForceKernel::CudaCalcShrinkForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, const System& system) :
        CalcShrinkForceKernel(name, platform), data(data) {
    for (int i = 0; i < (int) data.contexts.size(); i++)
        kernels.push_back(Kernel(new CommonCalcShrinkForceKernel(name, platform, *data.contexts[i], system)));
}

void CudaCalcShrinkForceKernel::initialize(const System& system, const ShrinkForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).initialize(system, force);
}

double CudaCalcShrinkForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    for (int i = 0; i < (int) data.contexts.size(); i++) {
        CudaContext& cu = *data.contexts[i];
        ComputeContext::WorkThread& thread = cu.getWorkThread();
        thread.addTask(new Task(context, getKernel(i), includeForces, includeEnergy, data.contextEnergy[i]));
    }
    return 0.0;
}

void CudaCalcShrinkForceKernel::copyParametersToContext(ContextImpl& context, const ShrinkForce& force) {
    for (int i = 0; i < (int) kernels.size(); i++)
        getKernel(i).copyParametersToContext(context, force);
}