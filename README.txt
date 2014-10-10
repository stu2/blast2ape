
    blast2ape v0.1
    Copyright 2014, Stuart Archer

    blast2ape is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    blast2ape is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can view the GNU General Public License at <http://www.gnu.org/licenses/>.




    Description:
    blast2ape is a simple script to import blast hits into ApE files so that they
    can be viewed as features in ApE. As this is meant to be used in conjuction
    with BLAST results from refseq RNA sequences, only sense-strand hits will be
    displayed.
    
    Be warned that ApE may run slowly or even crash if the minimum alignment size
    is too low, because thousands more hits will be displayed as features.
    
    Instructions:

    1) Make sure you have ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
       working on your system. blast2ape was built and tested with python 2.7.5.
    2) Place blast2ape.py, your downloaded blast hits table and your target sequence
       (as an .ape file) in the same directory.
    3) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
       'cd ' and then dragging the folder into the terminal from finder.
    4) Run the script by entering:

       python blast2ape.py  your_blast_output  your_ape_file  <minimum alignment length>  <maximum alignment length>  <min percent identity> <-i>'

       Items in <> are optional, but must be in order.

    5) Open the resultant file (blasthits_apefile) in ApE. Shorter hits (minimum
       alignment length) should be red, longer hits should be progressively orange
       then yellow. Remember, some long hits are hidden behind shorter hits, it is
       sometimes necessary to select a region and look at the features to see their
       actual size.

