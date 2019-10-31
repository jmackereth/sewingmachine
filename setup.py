from setuptools import setup

setup(name='sewingmachine',
      version='0.1',
      description='the Spectral Equivalent Widths(ing) machine',
      url='http://github.com/jmackereth/sewingmachine',
      author='J. Ted Mackereth',
      author_email='J.E.Mackereth@2011.ljmu.ac.uk',
      license='MIT',
      packages=['sewingmachine'],
      install_requires=['apogee', 'numpy', 'scipy', 'matplotlib', 'tqdm'],
      zip_safe=False)
