# VariantFormatter Operation Manual

## Operation
simpleVariantFormatter.format(<VARIANT(s)>, GENOMC_BUILD>, <TRANSCRIPT_SET>, <SPECIFY_TRANSCRIPT(s))

# Single variant processing
```
import json
import VariantFormatter
import VariantFormatter.simpleVariantFormatter
result_from_string = VariantFormatter.simpleVariantFormatter.format('15-72105928-AC-A', 'GRCh37', 'all', None)
print json.dumps(batch_results_from_string, sort_keys=False, indent=4, separators=(',', ': '))
```

# Batch processing
```
import json
import VariantFormatter
import VariantFormatter.simpleVariantFormatter
batch_results_from_string = VariantFormatter.simpleVariantFormatter.format('15-72105928-AC-A|NC_000012.11:g.122064777C>A|2-209113113-G-A,C,T', 'GRCh37', 'all', None)
print json.dumps(batch_results_from_string, sort_keys=False, indent=4, separators=(',', ': '))

Note - VF will also accept a variant list rather than a pipe delimited string
```

Accepted formats for variants are:
```
15-72105928-AC-A - dash delimited
15:72105928:AC:A - colon delimited
2-209113113-G-A,C,T - Multiple ALTS
NC_000012.11:g.122064777C>A - HGVS
```
Possible assemblies are:
```
GRCh37
GRCh38
```

You can specify a specific transcript set:
```
None or 'all' = RefSeq and Ensembl
refseq = RefSeq only
ensembl = Ensembl only
```

You can specify single or use multiple transcripts with: 
```
'NM_022356.3| NM_001146289.1| NM_001243246.1'
``

VariantFormatter produces a dictionary output that contain all possible HGVS descriptions
of the input variant, according to the user-specified parameters above 

## Unit testing

VariantFormatter is written to be pytest-compatible. Run
`pytest`
in the VariantFormatter root folder, the same as that in which this file resides. 
The test will take several minutes to complete, but runs through several hundred currated
variants

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