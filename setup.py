from setuptools import setup

setup(
    name='methylpy',
    version='0.1.1',
    author='Matthew Schultz',
    author_email='schultzmattd@gmail.com',
    packages=['methylpy', 'methylpy.test', 'methylpy.test.test_files'], 
    url='http://pypi.python.org/pypi/methylpy/',
    license='LICENSE.txt',
    description='A package for the processing and analysis of bisulfite sequencing data.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.6.1",
	"scipy >= 0.10.0",
    ],
    package_data = {
	    'methylpy.test':['test_files/*'],
	    'methylpy':['rms.cpp', 'run_rms_tests.out']
	}
)
