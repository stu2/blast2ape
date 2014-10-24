#!/usr/bin/env python

####################################################################################
#    xml2ape v0.1
#    Copyright 2014, Stuart Archer
#
#    xml2ape is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    xml2ape is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You can view the GNU General Public License at <http://www.gnu.org/licenses/>.
####################################################################################


####################################################################################
#
#    Instructions:
#
#
####################################################################################


import sys
try:
    import Bio
except:
    print "Cannot detect Biopython. Cannot calculate primer Tms\n"
    calc_tms = False
else:
    calc_tms =True
    import Bio
    from Bio.SeqUtils import MeltingTemp
import re
import os


usage = 'Usage: python ap.py  blast_out_xml_file  apefile  <minimum alignment length>  <maximum alignment length>  <min tm>'
   
apefile = sys.argv[1]
if not os.path.isfile(apefile):
    print apefile+' does not appear to exist in the current directory. Exiting.'
    print usage
    exit()

tablefile = sys.argv[2]

features = list()
query = str('')
with open(apefile, 'r') as ape_in:
    in_features=False
    in_sequence=False
    for line in ape_in:
        if line.startswith('FEATURES'):
            in_features = True
            continue
        if line.startswith('          '):
            continue
        if line.startswith('ORIGIN'):
            in_sequence = True
            in_features = False
            continue
        if in_features:
            mobj=re.match('\s+([\S]+)\s+([\d]+)\.\.([\d]+)$', line)
            if mobj:
                features.append([mobj.group(1), mobj.group(2), mobj.group(3)])
        elif in_sequence:
            line_seq = re.sub(r'[\d\s/]', '', line)  # take out numbers, whitespace
            query += line_seq.rstrip()
            
with open(tablefile, 'w') as tout:
    for row in features:
        if not row[0].startswith('blast_hit'):
            seq = (  query[int(row[1])-1:int(row[2])]  )
            row.append(seq)
            if calc_tms:
                tm = MeltingTemp.Tm_staluc(seq)
                row.append( tm )
            row = list(str(r) for r in row)
            print >> tout, '\t'.join(row)



