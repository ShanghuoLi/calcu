from setuptools import setup

setup(
  name = 'calcu',
  packages = ['calcu'],
  version = '0.1.0',
  description = 'To compute the physical parameters of molecular outflow and total column density of molecule',
  author = ['Shanghuo Li'],
  author_email = 'shanghuo.li@gmail.com',
  url = 'https://github.com/ShanghuoLi/calcu', # use the URL to the github repo
  keywords = ['astrophysics', 'calcu', 'outflow', 'column density'], # arbitrary keywords
  classifiers = [],
  install_requires=[
      'numpy',
      'matplotlib',
      'astropy',
  ]
)
