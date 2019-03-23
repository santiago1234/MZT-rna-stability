from setuptools import setup


def readme():
  with open('README.rst') as f:
    return f.read()


setup(name='mztrnastabilitty/',
      version='0.1',
      description='Python Code for the project 190108-mzt-rna-stability',
      long_description=readme(),
      url='http://github.com/storborg/funnies://github.com/santiago1234/MZT-rna-stability',
      author='Santiago Medina',
      author_email='smedina@stowers.org',
      license='MIT',
      packages=[
          'mztrnastabilitty',
          'mztrnastabilitty.dataprocessing'
      ],
      install_requires=[
          'pandas',
          'numpy',
          'argparse'
      ],
      scripts=[
          # does your project contain scripts
      ],
      zip_safe=False)
