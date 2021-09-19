from setuptools import Extension, setup

setup(
    name= 'mykmeanssp',
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['kmeans.c'],
        ),
    ]

)