from setuptools import setup, find_packages, Extension 

setup(
    name= 'mykmeanssp',
    install_requires=['invoke'],
    packages=find_packages(),
    ext_modules=[
        Extension(
            'spkmeans',
            ['spkmeans.c', 'spkmeansmodule.c'],
        ),
    ]
)

