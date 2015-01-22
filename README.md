Bioinformatics-Unitigging-Tool
==============================

Tool for unitigging created as part of coursework for Bioinformatics course on University of Zagreb, Faculty of Electrical Engineering and Computing (UNIZG FER) , academic year 2014/2015.

To compile and install program go to source directory and execute following.
```
make -f Makefile-1
```

To run program execute following command from command line.
```
./program <path-to-folder>
```

Folder where files with data are located has to contain following files:
* overlaps.afg - File with overlaps.
* reads.2k.10x.fasta - File which maps reads to it's name.
* reads.bnk/RED.0.map - File which maps every read name to it's id.
* 

All files are produced with Minimus and their description can be found on http://schatzlab.cshl.edu/teaching/2013/pacbioasm.shtml .

To test program go to test folder and execute
```
make -f Makefile
```

in command line.

Test program expects two arguments: path to layouts.afg file created by Minimus and layouts.afg file created by program.
Program compares layouts in both files and reports whether they are equal or not.
Example: 
```
./test ds/layouts.afg layouts.afg
```
