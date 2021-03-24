from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '/usr/local/openmm'
MSFplugin_header_dir = '/home/peanut/Software/molecular_shrink_ray/openmmapi/include'
MSFplugin_library_dir = '/home/peanut/Software/molecular_shrink_ray/cmake-build-debug'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

openmm_lib_path = os.getenv('OPENMM_LIB_PATH')

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_lib_path]

extension = Extension(name='_MSFPlugin',
                      sources=['MSFPluginWrapper.cpp'],
                      libraries=['OpenMM', 'MSFPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), MSFplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), MSFplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                      )

setup(name='MSFPlugin',
      version='1.0',
      py_modules=['MSFPlugin'],
      ext_modules=[extension],
      )
