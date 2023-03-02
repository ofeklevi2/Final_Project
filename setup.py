from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['spkmeansmodule.c', "spkmeans.c"], extra_compile_args = ["-O0"])
setup(name='mykmeanssp',
     version='1.0',
     description='Python wrapper for spkmeans C extension',
     ext_modules=[module])