# -*- coding:utf-8;mode:python -*-
# Copyright [2019] tadashi9e
from distutils.core import setup, Extension

setup(name='qc', version='0.0.2',
      packages=['qc'],
      ext_modules=[Extension('qc_py', ['qc_py.cc', 'qc.cc'],
                             extra_compile_args=['-std=c++11'])])
