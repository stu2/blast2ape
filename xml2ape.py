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
#    1) Make sure you have ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
#       working on your system. xml2ape was built and tested with python 2.7.5.
#    2) Place xml2ape.py, your downloaded blast hits table and your target sequence
#       (as an .ape file) in the same directory.
#    3) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
#       'cd ' and then dragging the folder into the terminal from finder.
#    4) Run the script by entering:
#
#       python xml2ape.py  your_blast_output  your_ape_file  <minimum alignment length>  <maximum alignment length>  <min percent identity> <-i>'
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
import Bio
import re
import numpy
import os
from Bio.SeqUtils import MeltingTemp
from Bio.Blast import NCBIXML
from Bio.Blast import Record

usage = 'Usage: python xml2ape.py  blast_out_xml_file  apefile  <minimum alignment length>  <maximum alignment length>  <min tm>'

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
    print 'Warning: fourth argument does not appear to be an integer. Assigning default maximum \n length (1000 nt or query length minus 2, whichever is shorter).'
    maxlen = 1000
else:
    maxlen = int(sys.argv[4])
    
try:
    mintm = float(sys.argv[5])
except:
    print 'Warning: fourth argument is absent or not a number. Assigning default minimum percent identity (100).'
    mintm = float(40.0)
else:
    mintm = float(sys.argv[5])
    
verbose = False
colour_scale = 15


from Bio import SeqIO
for seq_record in SeqIO.parse(apefile, "genbank"):
    full_query = str(seq_record.seq)
### print query sequence: ###
#    print full_query
    qlen = len(full_query)
    if qlen - 2 < maxlen:
        maxlen = qlen - 2
        print "warning: BLAST query was shorter than maximum hit length. To avoid identifying its"
        print "own transcript as an off-target, setting maxlen to 2 nt less than the query length."

opp = int(0)
same = int(0)
cross_hits=dict()
with open(bfile, 'r') as blin:
    brecs = NCBIXML.read(blin)
    for alignment in brecs.alignments:
        for hsp in alignment.hsps:

            # skip if opposite strand:
            if hsp.sbjct_end < hsp.sbjct_start:
                opp += 1
                continue
            same += 1
            
            query = hsp.query
            match = hsp.match
            transcript = alignment.title   
            spaced = ""
            for i in range(len(match)):
                if match[i] == " ":
                    spaced += "-"
                else:
                    spaced += hsp.query[i]
            for seq in spaced.split('-'):
                if len(seq) < minlen:
                    continue
                if len(seq) > maxlen:
                    continue
                tm = MeltingTemp.Tm_staluc(seq)
                if tm < mintm:
                    continue
                if seq in cross_hits:
                    print "warning: more than one blast hit with identical cross-homology sequence. Reporting one only."
                    cross_hits[seq][2] += 1
                    continue
                cross_hits[seq]=[tm, transcript, 1]
print "same: "+str(same) + " opp: "+str(opp)


seqs = list(seq for seq in cross_hits)  # get unordered list of seqs
tms = numpy.array(list(cross_hits[seq][0] for seq in seqs)) # get list of tms in same order as the list of seqs
temp = tms.argsort()
ranks = numpy.empty(len(tms), int)
ranks[temp] = numpy.arange(len(tms)) # this is the rank for both tms and seqs (each rank is unique even if tm is a tie)
ordered_seqs = ['na'] * len(ranks)
for i in enumerate(ranks, 0):
    ordered_seqs[i[1]] = seqs[i[0]]

feats = str()
for seq in ordered_seqs:  # make high-tm hits last; thus they are displayed on top of the others
    tm, transcript, multi = cross_hits[seq]
    if verbose:
        print "tm: "+str(tm)+"number of off-target transcripts with this hit: "+str(multi)+ "reported transcript: "+transcript
    green = format(int(min((tm-mintm)*colour_scale, 255)), '02x')
    colour = '#ff'+str(green[:2])+'00'
    starts = [m.start() for m in re.finditer('(?=%s)' % seq, full_query)]
    for st in starts:
        end = st + len(seq)
        feats += '     blast_hit       ' +str(st+1)+'..'+str(end)
        feats += '\n                     /label=tm%.2f' % tm  +"|"+transcript[:42]
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
print 'xml2ape.py output in: blasthits_'+apefile+'. \n'


