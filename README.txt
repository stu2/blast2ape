
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
    
    blast2ape contains three python scripts:
    xml2ape.py (for incorporating BLAST output into ApE files)
    ape2table.py (for exporting probes designed in ApE)
    blast2ape.py (a simpler alternative to xml2ape.py)
    
    
    
    #########################################################
    #                   xml2ape.py                          #
    #########################################################    

    Description:
    xml2ape.py is a script for parsing BLAST output and assisting its display in the
    graphical plasmid editing software, ApE. Regions of cross-homology, identified
    through a BLAST search and downloaded as an XML file, can be displayed.
    
    Dependencies:
    xml2ape.py requires:
    * Biopython (www.biopython.org). Blast2ape.py (below) is an alternative script to
    use if Biopython is unavailable, however it is not recommended.
    * Python2.7.5
    * ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
    
    Instructions:
    
    1) Make sure you have ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
       working on your system.
    2) Place xml2ape.py, your downloaded blast hits xml file and your target sequence
       (as an .ape file) in the same directory.
    3) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
       'cd ' and then dragging the folder into the terminal from finder.
    4) Run the script by entering:

    python xml2ape.py -b blast_out_xml_file -a apefile  <-l minimum alignment length>
    <-x maximum alignment length>  <-t min tm> <-r>'
    
    
    -b : blast output xml file, downloaded after BLAST search.
    -a : A pre-existing .ape file containing the target sequence (the query used in the
        BLAST search.)
    -l : minimum alignment length. Default: 15
    -x : Maximum alignment length. Without a maximum hit length, xml2ape would recognize
        the target transcript itself as an off-target hit and display it. Default: 1000
    -r : Binary flag. Include reverse-strand hits (not recommended for making pdd probes)
    -t : minimum tm of hits to report. Calculated as by Allawi and SantaLucia, 1997
        using the default parameters as implemented in Biopython.
    Items in <> are optional.


    #########################################################
    #                   ape2table.py                        #
    #########################################################    
    
    Description:
    ape2table is a python script for converting features that have been designed in APE to table format. 
    
    Instructions:
    
    1) Place ape2table.py and your .ape file in the same directory.
    2) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
       'cd ' and then dragging the folder into the terminal from finder.
    3) Run the script by entering:
    
    python ape2table.py  -a apefile.ape <-o output_file.txt> <-i>
    
    explanation of options:

    -a  input apefile containing features to extract
    -o  output table filename. Output will be tab-delimited text.
    -i  option not to filter out features beginning with "blast_hit"
    Items in <> are optional.
    

    #########################################################
    #                   blast2ape.py                        #
    #########################################################
    

    Description:
    blast2ape.py is intended as an alternative to xml2ape.py if Biopython is not
    installed on the system, however the output is not equivalent.
    
    blast2ape.py is a simple script to import blast hits into ApE files so that they
    can be viewed as features in ApE. As this is meant to be used in conjuction
    with BLAST results from refseq RNA sequences, only sense-strand hits will be
    displayed. If you have Biopython installed, it is preferable to use xml2ape.py
    (an accompanying script, see above) to parse BLAST results as it splits up
    hits into stretches of perfect homology for display, and also calculates match
    tm, while blast2ape displays these hits as one long feature provided they exceed
    minimum % identity thresholds, and does not calculate tm.
    
    Be warned that ApE may run slowly or even crash if the minimum alignment size
    is too low, because thousands more hits will be displayed as features. Also note
    that the number of hits reported may be more than in the XML file, as imperfectly
    matched hits will be split up at the site of the mismatches into smaller hits.
    
    
    Instructions:
    
    ** Note that blast2ape.py takes BLAST output that has been downloaded as a table,
    not an XML file **

    1) Make sure you have ApE (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
       working on your system. blast2ape was built and tested with python 2.7.5.
    2) Place blast2ape.py, your downloaded blast hits table and your target sequence
       (as an .ape file) in the same directory.
    3) Open up a terminal, cd to the directory. In Mac OSX, you can do this by typing
       'cd ' and then dragging the folder into the terminal from finder.
    4) Run the script by entering:

       python blast2ape.py  your_blast_output  your_ape_file  <minimum alignment length>
       <maximum alignment length>  <min percent identity> <-i>'

    Items in <> are optional, but must be in order. if '-i' is added, the script will
    ignore the field order stated in the header
       

    5) Open the resultant file (blasthits_apefile) in ApE. Shorter hits (minimum
       alignment length) should be red, longer hits should be progressively orange
       then yellow. Remember, some long hits are hidden behind shorter hits, it is
       sometimes necessary to select a region and look at the features to see their
       actual size.
       

