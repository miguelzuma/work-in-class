# work-in-class
Suite of utilities to work with the CLASS and HI_CLASS codes

Note it is under heavy development and may contain bugs and/or be useless for you (then fork, write and pull-request!)

You need to follow one of these procedures in order to use this package: 

1. Automate installation (recommended):

- From github:

        pip install --user git+https://github.com/ardok-m/work-in-class

    or

        pip install --user git+https://github.com/miguelzuma/work-in-class

(normally ardok-m has more recent features but also is more prone to bugs)

- After downloading the package:

        pip install --user /path/to/package

    or

        python /path/to/package/setup.py install --user

2. Source the folder where the utilities are located:

- for ipython (notebook) add

        import sys
        sys.path.insert(0, '/path/to/package/work_in_class')

- for python (e.g. in terminal) add to .bashrc

        export PYTHONPATH=$HOME/path/to/package/work_in_class/:$PYTHONPATH

To Do (Sept 22 '16):
- Add functions to compare output, generate ini files
- Add/test compatibility with Classy
- Load and handle init files. Add a dictionary for abbreviations to parameters and nice prints.
- Read the headers and automatically translate column numbers into labels

** This list is unmantained at the time being.
