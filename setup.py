
from setuptools import setup, find_packages
__version__='1.0.1'
__package__='coxmunk'

setup(
    name=__package__,
    version=__version__,
    packages=find_packages(exclude=['build']),
    package_data={'':['*.so','*h','*angles*']},
    #     # If any package contains *.txt files, include them:
    #     '': ['*.txt'],
    #     'lut': ['data/lut/*.nc'],
    #     'aux': ['data/aux/*']
    # },
    include_package_data=True,

    url='',
    license='MIT',
    author='T. Harmel',
    author_email='tristan.harmel@gmail.com',
    description='simple sunglint computation based on Cox-Munk wave slope distribution',

    # Dependent packages (distributions)
    install_requires=['pandas','numpy','netCDF4','xarray',
                      'matplotlib','docopt'],

    entry_points={
          'console_scripts': [
              'coxmunk = coxmunk.visu:main'
          ]}
)
