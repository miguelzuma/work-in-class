import os
from distutils.core import setup

# TODO: consider using setuptools instead of distutils:
# https://packaging.python.org/guides/tool-recommendations/


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="work-in-class",
    version="dev",
    description=("Tools to work with hi_class and MontePython"),
    license="GPL3",
    keywords="class hi_class montepython",
    url="https://github.com/miguelzuma/work-in-class",
    packages=['work_in_class', 'work_in_class.plot_class'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics"
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
