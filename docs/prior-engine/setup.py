#!/usr/bin/env python

from setuptools import setup

with open('requirements.txt') as r:
    requirements = r.read().splitlines()
requirements = [r.split("=")[0] for r in requirements]

__version__ = None
__status__ = None
__license__ = None

with open('multiply_prior_engine/version.py') as f:
    exec(f.read())

setup(name='multiply-prior-engine',
      version=__version__,
      description='MULTIPLY Prior Engine',
      author='MULTIPLY Team',
      packages=['multiply_prior_engine'],
      package_data={
          'multiply_prior_engine': ['*.yml']
      },
      entry_points={
          'prior_creators': [
              'vegetation_prior_creator = multiply_prior_engine:vegetation_prior_creator.VegetationPriorCreator',
              'soil_moisture_prior_creator = multiply_prior_engine:soilmoisture_prior_creator.SoilMoisturePriorCreator',
              'roughness_prior_creator = multiply_prior_engine:soilmoisture_prior_creator.RoughnessPriorCreator',
          ],
          'console_scripts': [
              'user_prior = multiply_prior_engine.user_prior:main'
          ]
      },
      install_requires=requirements
      )
