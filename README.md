# coxmunk
Simple sunglint computation based on Cox and Munk model





## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

What things you need to install the software and how to install them

```
python3 -m pip install --user --upgrade setuptools
```

### Installing

First, clone [the repository](https://gitlab.irstea.fr/telquel-obs2co/insitu/trios#) and execute the following command in the
local copy:

```
python3 setup.py install 
```

This will install the package into the system's Python path.
If you have not the administrator rights, you can install the package as follows:

```
python3 setup.py install --user
```

If another destination directory is preferred, it can be set by

```
python3 setup.py install --prefix=<where-to-install>
```

This installation is supposed to download
and compile all the associated packages as well as prepare the executables `trios_processing` and `trios_visual`.

If the installation is successful, you should have:
```
$ trios_processing
Usage:
  trios_processing <input_dir> <IDpr> <measurement_type> --lat <lat> --lon <lon>    [--altitude=alt] [--ofile <ofile>] [--odir <odir>] [--plot] [--figdir <figdir>]    [--name <name>] [--method <method>] [--no_clobber]
  trios_processing -h | --help
  trios_processing -v | --version
```

## Running the tests
