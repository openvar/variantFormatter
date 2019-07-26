# VariantFormatter Installation

You will need to install [VariantValidator](https://github.com/openvar/variantValidator) >1.0 and the necessary databases identified in the INSTALLATION.md

The packages required for VariantValidator to function are now set up in the environment "vvenv".

## Installing VariantFormatter

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
