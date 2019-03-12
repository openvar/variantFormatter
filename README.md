                         # About VariantFormatter

VariantValidator is a software tool designed to accept genomic level sequence variant 
descriptions in the pseudo-VCF or HGVS formats (see MANUAL.md), and automatically generate 
relevant transcript-level and protein-level HGVS variant descriptions.

VariantFormatter interfaces with a custom version of the hgvs package to parse, format, 
and manipulate biological sequence variants.  
See https://github.com/biocommons/hgvs/ for details of the hgvs package
See https://github.com/openvar/vv_hgvs for the custom version of the hgvs package


## Features

The basic functionality of VariantFormatter is documentes in the MANUAL.md

VariantFormatter simultaneously and accurately projects genomic sequence variations onto 
all overlapping transcript reference sequences

VariantFormatter can be used as a variant-by-variant processor or as a batch processor

Alternatively, genomic sequence variation can be projected onto a specified single, or 
specified subset of transcript reference sequences for any given gene

Projection of sequence variations between reference sequences takes account of 
discrepancies between genomic and transcript reference sequences, thus ensuring an 
accurate prediction of the effect on encoded proteins for every gene

## Web-API
http://rest.variantvalidator.org

## Pre-requisites

VariantFormatter will work on Mac OS X or Linux-compatiable computers.

Required software:
* MySQL
* Python 2.7

Optional software:
* Postgres version 9.5 or above
* SQLite version 3.8.0 or above

For installation instructions please see [INSTALLATION.md](INSTALLATION.md)

# Operation and configuration

Please see [MANUAL.md](MANUAL.md)

## License

Please see [LICENSE.txt](LICENSE.txt)

## Cite us
VariantFormatter is a sub-project of the VariantValidator project

Hum Mutat. 2017 Oct 1. doi: 10.1002/humu.23348

VariantValidator: Accurate validation, mapping and formatting of sequence variation 
descriptions

Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R.

> Copyright (C) 2019 Peter Causey-Freeman, University of Leicester
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
