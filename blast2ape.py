#!/usr/bin/env python

####################################################################################
#    blast2ape v0.1
#    Copyright 2014, Stuart Archer
#
#    blast2ape is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    blast2ape is distributed in the hope that it will be useful,
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
#    1) Make sure you have ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
#       working on your system. blast2ape was built and tested with python 2.7.5.
#    2) Place blast2ape.py, your downloaded blast hits table and your target sequence
#       (as an .ape file) in the same directory.
#    3) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
#       'cd ' and then dragging the folder into the terminal from finder.
#    4) Run the script by entering:
#
#       python blast2ape.py  your_blast_output  your_ape_file  <minimum alignment length>  <maximum alignment length>  <min percent identity> <-i>'
#
#       Items in <> are optional, but must be in order.
#
#    5) Open the resultant file (blasthits_apefile) in ApE. Shorter hits (minimum
#       alignment length) should be red, longer hits should be progressively orange
#       then yellow. Remember, some long hits are hidden behind shorter hits, it is
#       sometimes necessary to select a region and look at the features to see their
#       actual size.
#
####################################################################################


import sys
import os
import math
usage = 'Usage: python blast2ape.py  blast_out_file  apefile  <minimum alignment length>  <maximum alignment length>  <-i>'

bfile = sys.argv[1]
if not os.path.isfile(bfile):
    print bfile+' does not appear to exist in the current directory. Exiting.'
    print usage
    exit()
    
apefile = sys.argv[2]
if not os.path.isfile(apefile):
    print apefile+' does not appear to exist in the current directory. Exiting.'
    print usage
    exit()
    
try:
    minlen = int(sys.argv[3])
except:
    print 'Warning: third argument does not appear to be an integer. Assigning default minimum length (15 nt).'
    minlen = 15
else:
    minlen = int(sys.argv[3])
    
try:
    maxlen = int(sys.argv[4])
except:
    print 'Warning: fourth argument does not appear to be an integer. Assigning default maximum length (1000 nt).'
    maxlen = 1000
else:
    maxlen = int(sys.argv[4])
    
try:
    min_pci = float(sys.argv[5])
except:
    print 'Warning: fourth argument is absent or not a number. Assigning default minimum percent identity (100).'
    min_pci = float(100.0)
else:
    min_pci = float(sys.argv[5])
    
default_fns = 'query id, subject ids, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score'   
try:
    str(sys.argv[6])
except:
    print "No fifth arg given, using field order defined in BLAST results header."
    header_defined_fields = 1
else:
    if str(sys.argv[6]) == '-i':
        header_defined_fields = 0
        print "Ignoring field order defined in BLAST results header, using default order:"
        print default_fns
    else:
        print "Fifth argument not \'-i\', using field order defined in BLAST results header."
    

feats = str()
i = 0
with open( bfile, 'r') as blin:
    for line in blin:
        if line.startswith('# Fields: '):
            fieldnames = line.rstrip('\n').rstrip('\r').split(', ')
            fieldnames[0] = fieldnames[0][10:]
            if not ', '.join(fieldnames) == default_fns:
                print 'Warning: fields are in an unexpected order in this BLAST results file.'
                print 'This version: '+', '.join(fieldnames)
                print 'Normal order: '+default_fns
                print 'Warning: This script will generate key errors and fail if the required fields are absent in the header.'
                print 'To force the expectation of fields in the normal order, ignoring this header line, input \'-i\' as the fifth argument to the script.'
        if line.startswith('#'):
            continue
        try:
            fieldnames
        except:
            print 'No field-names given in header. Using defaults'
            fieldnames = default_fns.split(', ')
        try:
            idx
        except:
            idx = dict()
            for f in fieldnames:
                idx[f] = fieldnames.index(f)
        splut = line.rstrip('\n').rstrip('\r').split('\t')
        if len(splut) < len(idx):
            continue
        # skip if opposite strand
        if int(splut[idx['s. start']]) > int(splut[idx['s. end']]):
            continue
        al_len = int(splut[idx['alignment length']])
        if al_len < minlen:
            continue
        if al_len > maxlen:
            continue
        if float(splut[idx['% identity']]) < min_pci:
            continue
        len_diff_log2 = math.log(1 + al_len - minlen, 2)
        green = format(int(min(len_diff_log2*100, 255)), '02x')
        colour = '#ff'+str(green[:2])+'00'
        
        feats += '     blast_hit       '+splut[idx['q. start']]+'..'+splut[idx['q. end']]
        feats += '\n                     /label='+splut[idx['subject ids']]+splut[idx['q. start']]+'..'+splut[idx['q. end']]
        feats += '\n                     /ApEinfo_fwdcolor='+colour
        feats += '\n                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}'
        feats += '\n'
        
with open(apefile, 'r') as ape_in:
    with open('blasthits_'+apefile, 'w') as ape_out:
        features_exist = 0;
        for line in ape_in:
            if line.startswith('FEATURES             Location/Qualifiers'):
                features_exist = 1
            if line.startswith('ORIGIN'):
                if not features_exist:
                    print >> ape_out, 'FEATURES             Location/Qualifiers'
                print >> ape_out, feats+'\n'
            print >> ape_out, line
print 'blast2ape.py output in: blasthits_'+apefile+'. \n'
