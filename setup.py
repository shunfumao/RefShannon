# python setup.py install --user

from setuptools import setup, find_packages
setup(
  name = 'RefShannon',
  version = '0.0.1',
  description = 'RefShannon',
  url='https://github.com/shunfumao/RefShannon',
  author='Shunfu Mao',
  install_requires=['pandas'
                   ],
  packages = find_packages()
)