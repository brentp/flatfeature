from setuptools import setup

setup(
    name = "flatfeature",
    version = "0.1.0",
    author = "Brent Pedersen",
    author_email = "bpederse@gmail.com",
    description = ("simple format for genomic features"),
    license = "BSD",
    keywords = "bioinformatics genomics bio",
    url = "http://packages.python.org/flatfeature",
    install_requires = ['numpy>=1.4.0', 'pyfasta'],
    packages=['.'],
    long_description=open('README.rst').read(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
