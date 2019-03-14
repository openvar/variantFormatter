# VarianFormatter Installation

These instructions will allow you to configure the software on Linux and Mac OS X computers.

There are several steps involved in setting up VariantFormatter:
* The python environment must be set up with the correct packages
* The VariantFormatter files themselves must be downloaded and installed.
* The databases must be downloaded and set up
* The configuration files must be changed to point VariantFormatter at those databases.

## Virtual environment (Python 2.7)

When installing Variant Validator it is wise to use a virtual environment, as it requires specific versions of several libraries.
We recommend using conda.
```
$ conda create -n VFenv
$ conda activate VFenv
$ conda install -c conda-forge sqlite python=2.7 protobuf=3.5.1
```
The packages required for variant validator to function are now set up in the environment "VVenv".

## Installing validator code

To clone this software from GIT, use:
```
$ git clone https://github.com/openvar/variantFormatter.git
$ cd variantFormatter/
```

Run the installation script to integrate VariantFormatter with python's site packages.
```
$ python setup.py install
```

For development purposes, you can use
```
$ pip install -e .
```
to ensure any changes you make in the local variant validator folder is reflected in your python site-packages.


## Setting up UTA database (PostGreSQL >=9.5)

It's recommended for performance reasons to use a local version of the UTA database. We again recommend creating a specific user account.
```
CREATE ROLE uta_admin WITH CREATEDB;
ALTER ROLE uta_admin WITH LOGIN;
\password
CREATE DATABASE uta WITH OWNER=uta_admin TEMPLATE=template0;
```

To fill this database, download the gzipped uta genetics database, and upload it into psql.
```
$ wget http://dl.biocommons.org/uta/uta_20180821.pgd.gz
$ gzip -cdq uta_20180821.pgd.gz | psql -U uta_admin -v ON_ERROR_STOP=0 -d uta -Eae
```


## Setting up Seqrepo (SQLite >=3.8)

VariantFormatter requires a local SeqRepo database. The seqrepo library is already installed, but you'll need to download a seqrepo database. These instructions assume you are using your home directory; however, you can place the database it so long as you modify the config file which VariantFormatter places in your home directory when prompted.
```
$ mkdir seqrepo
$ seqrepo --root-directory ~/seqrepo pull -i 2018-08-21
```
To check it has downloaded:
```
$ seqrepo --root-directory ~/seqrepo list-local-instances
```

## Configuration

See the [manual](MANUAL.md) for configuration instructions.


> <LICENSE>
> Copyright (C) 2019  Peter Causey-Freeman, University of Leicester
> 
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU Affero General Public License as
> published by the Free Software Foundation, either version 3 of the
> License, or (at your option) any later version.
> 
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU Affero General Public License for more details.
> 
> You should have received a copy of the GNU Affero General Public License
> along with this program.  If not, see <https://www.gnu.org/licenses/>.
> </LICENSE>
