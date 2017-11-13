from setuptools import setup

setup(
    name='methylpy',
    version='1.1.1',
    author='Yupeng He',
    author_email='yupeng.he.bioinfo@gmail.com',
    packages=['methylpy'], 
    url='http://pypi.python.org/pypi/methylpy/',
    license='LICENSE.txt',
    description='Bisulfite sequencing data processing and differential methylation analysis',
    long_description=open('docs/README.rst').read(),
    install_requires=[
        "numpy >= 1.6.1",
	"scipy >= 0.10.0",
    ],
    package_data = {
        'methylpy.test':['test/*'],        
	'methylpy':['rms.cpp', 'run_rms_tests.out']
    },
    keywords = ["Bioinformatics pipeline",
                "DNA methylation",
                "Bisulfite sequencing data",
                "Nome-seq data",
                "Differential methylation",
                "Calling DMRs",
                "Epigenetics",
                "Functional genomics"],
    scripts = ["bin/methylpy"]
)
