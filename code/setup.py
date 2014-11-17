#!/usr/bin/env python2

from setuptools import setup, Extension
import os

libraries = [
    'gmp',
    'gomp',
]
compile_args = [
    '-fopenmp',
    '-O3',
    '-Wall',
    '--std=gnu99',
    '-Isrc/gghlite-flint',
    '-Isrc/gghlite-flint/dgs',
    '-Isrc/gghlite-flint/dgs/dgs',
    '-Isrc/gghlite-flint/flint',
    '-DFLINT_CPIMPORT="src/gghlite-flint/flint/qadic/CPimport.txt"',
]

gghlite_srcs = []
for root, dirs, files in os.walk('src/gghlite-flint'):
    for file in files:
        if file.endswith('.c'):
            gghlite_srcs.append(os.path.join(root, file))
gghlite_srcs = filter(lambda s: '/benchmarks/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/doc/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/examples/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/link/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/profile/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/test/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/tests/' not in s, gghlite_srcs)
gghlite_srcs = filter(lambda s: '/applications/' not in s, gghlite_srcs)

zobfuscator = Extension(
    'obf._zobfuscator',
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/circuit.c',
        'src/clt_mlm.c',
        'src/_zobfuscator.c',
        'src/mpn_pylong.c',
        'src/mpz_pylong.c',
        'src/pyutils.c',
        'src/utils.c',
    ]
)

obfuscator = Extension(
    'obf._obfuscator',
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/clt_mlm.c',
        'src/clt_utils.c',
        'src/ggh_mlm.c',
        'src/_obfuscator.c',
        'src/mpn_pylong.c',
        'src/mpz_pylong.c',
        'src/pyutils.c',
        'src/utils.c',
    ] + gghlite_srcs
)

setup(name='obfuscator',
      author='Alex J. Malozemoff',
      author_email='amaloz@cs.umd.edu',
      version='1.0alpha',
      description='Implementation of cryptographic program obfuscation',
      license='GPLv2',
      url='https://amaloz.github.io/obfuscation',
      packages=['obf'],
      ext_modules=[obfuscator, zobfuscator],
      scripts=['obfuscator'],
      test_suite='t',
      classifiers=[
          'Topic :: Security :: Cryptography',
          'Environment :: Console',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Operating System :: Unix',
          'License :: OSI Approved :: Free For Educational Use',
          'Programming Language :: C',
          'Programming Language :: Sage',
      ])
