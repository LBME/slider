# SLIDER_VENOM.py

20 September 2021
Author : Rafael Borges
Objective: Generate amino acid possibilities by residue, fit best rotamer in electron density and calculate their
side chain and main chain real-space correlation coefficient

```
USAGE: SLIDER_VENOM.py pdb_file reflections.mtz map_coeffs.mtz output TRYALL
```

pdb_file should be a protein coordinate file with best possible fit against electron density
reflections.mtz should contain intensities or amplitudes
map_coeffs.mtz with calculated map (2FOFCWT/PH2FOFCWT)
output should be a new path, a folder 'output' will contain all SLIDER output and files starting with output will
    summarize information:
type_of_run (which was TRYALL in example USAGE)

Type_of_run should be a string with keywords and can change how SLIDER is run:
Amino acid possibilities can be generated in three different scenarios:
Trying all 20 possibilities for each residue (keyword: TRYALL)
Possibilities restricted by:
 mass spectrometry (keyword: MASSPEC)
 alignment (keyword: ALIGN)
If the last two options are chosen, an additional file containing either mass spectrometry or alignment file should be given.

The mass spectrometry file should be a file containing text of amino acids (aa) in one letter code (FASTA), each column should contain generated aa, per example, if 1st and 3rd residue should be a F and T, and 2nd residue either L,D or N, the file should be:

```
FLT
D
N
```

Keyword SKIPTEST should be given to skip RAM memory calculation;
 
Links to the external software of Phenix and Coot should be accessible through the terminal.
 
Output files:
output_all.log contains information of chain and residue number, RSCC and delta contrast
output_Resolved.log same as output_all.log, but only resolved residues
output_Dubious.log  same as output_all.log, but only dubious residues
output_runline.log has the line used to run SLIDER_VENOM.py
output_polder_coot_open_maps contains a script to open omit maps in coot (Go to Calculate -> Run Script -> select file)





# ConSurf_query_msa_SLIDER.py

19 April 2020
Author : Rafael Borges
Objective: Given ConSurf query_msa.aln and SLIDER table,
return percentage of SLIDER residue

```
USAGE: ConSurf_SLIDER.py query_msa.aln output_Resolved.log OutputFile
```

OutputFile should be a new path

RETURNS: a table with SLIDER details plus absolute number and % of 100% of given residue




#GivenDubiousResolvedResiduesReturnSequenceVariability.py

20 September 2021
Author : Rafael Borges
Objective: Given Dubious and Resolved log files from SLIDER_VENOM with residues and their RSCC,
return variability of sequences

```
USAGE: GivenDubiousResolvedResiduesReturnSequenceVariability.py Dubious.log Resolved.log OutputFile
```

OutputFile should be a new path

RETURNS: a file with amino acid possibilities (line) by residue number (column)



#ObtainPeptidesPatternLabRun.py

08 April 2020
Author : Rafael Borges
Objetive: Obtain positive SLIDER peptides from overall run from PATTERNLAB
It reads a sequence and patternlab fasta file with all values
It returns each peptide among its PrimaryScore and PPM

```
USAGE: ObtainPeptidesPatternLabRun.py seq.seq fasta.fasta OutputFile
```

seq.seq has just the sequence without >
fasta.fasta exactly like PatternLab output
OutputFile should be a new path

RETURNS: a table with peptide PrimaryScore PPM
         a file with each peptide with its position in sequence





#ObtainPeptidesPEAKSRunHTML.py

08 April 2020\
Author : Rafael Borges\
Objetive: Obtain positive SLIDER peptides from overall run from PEAKS HTML output\
It reads a sequence and PEAKS HTML output with all values\
It returns each peptide among its -10lgP\
```
USAGE: ObtainPeptidesPatternLabRun.py seq.seq PEAKS.html OutputFile
```
seq.seq has just the sequence without >\
PEAKS.html exactly like PEAKS html output\
OutputFile should be a new path\
\
RETURNS: a table with peptide PrimaryScore PPM\
         a file with each peptide with its position in sequence\
\
(C) 2021 Rafael Junqueira Borges 
