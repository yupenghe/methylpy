from setuptools import setup

methylpy_version = '1.4.7'

setup(
    name='methylpy',
    version=methylpy_version,
    author='Yupeng He',
    author_email='yupeng.he.bioinfo@gmail.com',
    packages=['methylpy'], 
    license='LICENSE.txt',
    description='Bisulfite sequencing data processing and differential methylation analysis',
    long_description=open('docs/README.rst').read(),
    install_requires=[
        "numpy >= 1.6.1",
	"scipy >= 0.10.0",
        "pysam >= 0.5.3"
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

if __name__ == '__main__':
    f = open("methylpy/__init__.py",'w')
    f.write("__version__ = \'"+methylpy_version+"\'"+"\n")
    f.close()
