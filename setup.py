from setuptools import setup, find_packages

setup(
  name = 'calcu',
  version = '0.1.0',
  description = 'To compute the physical parameters of molecular outflow and total column density of molecule',
  author = ['Shanghuo Li'],
  author_email = 'shanghuo.li@gmail.com',
  url = 'https://github.com/ShanghuoLi/calcu', # use the URL to the github repo  
  packages=find_packages(),
  classifiers = [
      'numpy',
      'matplotlib',
      'astropy',
  ]
)
