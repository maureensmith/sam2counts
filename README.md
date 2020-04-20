# sam2counts

The initial purpose of sam2counts was to generate count files from resulting MIME sequencing data which are used as input for
the [MIMEAnTo](https://github.com/maureensmith/MIMEAnTo) tool. 
It can be used in general to count nucleotides occurrences for each position (-pair or -triplet) with paired-end read
mappings in SAM format. 



The program is called from the command line with
```
sam2countsProg -r <reffile> -s1 <samfile1> [-s2 <samfile2>] -o <outfile> -d <dimension> [-q <qualiThreshold>] [-a]
```
and requires the following parameters

| flag      | parameter       | type          | description  |
| :---------| :-------------  |:-------------:| :-----|
|-r, -ref      | reffile         | (string)      |   path to reference sequence file in fasta format |
|-s1, -sam1    | samfile1        | (string)      |   path to (first) SAM file of the mapped reads |
|-s2, -sam2    | samfile2        | (string)      |   (optional) path to second SAM file of the mapped paired reads |
|-o, -out      | outfile         | (string)      |   path to output file containing the counts in tsv format |
|-d, -dimension| dimension       | (int)         |   either 1, 2, or 3 for the counting dimension  |
|-q, -quality  | qualiThreshold  | (int)         |   (optional) threshold for the quality score |
|-a, -ambig    |                 |               |   (optional) flag if ambiguous symbols shall be counted|

By default, only the nucleotides are counted and any ambiguity is ignored. If this flag ist set, 
also all ambiguity symbols according to the IUPAC notation are counted. For now, this option is 
only enabled for dimension 1 and for single read SAM files. 

We added a quality threshold despite the option of quality clipping in preprocessing steps. 
The reason is, that it is usually assumed that the quality of a
sequence drops towards the ends of the fragments. Hence, nucleotides with bad quality
are trimmed until a good quality is reached. However, it may also happen that a drop
of quality occurs within the sequence and these nucleotides of bad quality remain in the
read. Therefore, we check for the quality of all nucleotides, since this requires only one
additional operation per iteration.
If a '*' is given in the quality field in the SAM or no quality threshold is given, the bases are counted without any quality check.

The program can be divided into three procedures. The read pairs are imported and prefiltered,
the cohesive sequences of a read pair are merged to one read with certain rules,
and finally, the reads are counted for each position and nucleotide and written into the
output file. The procedures are described in particular below. 

## Prefilter
After importing the reference sequence, the read pairs are imported from the two SAM
files and processed one by one. Each alignment line of a SAM file contains 11 mandatory
fields of essential mapping information, and a variable number of optional fields. An
example of an alignment line is given below:
```
FCC1073ACXX:6:1101:17307:2200#AGTCAAAT/1 0 HIV1_535 22 70 94M6S * 0 0
TCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTTGTGGT
abbccccegggggiiiiihhhhiiiihhiihhhgfghghfffiiiffhiiiifihicfhffhhiiiheggdgeeeceeadb_bddccbbbaaaBBBBBBB
PG:Z:novoalign AS:i:64 UQ:i:64 NM:i:0 MD:Z:94
```
Only those fields, which are used for further filtering are explained here.
In a first filtering step, both sequences are checked for their mapping status. We considered
multiple conditions to assess mapped reads:
* The second field in a line contains a bit flag. If the fourth bit is set, the fragment is
unmapped.
* The fourth field gives the left most mapping position. To be mapped, the positions
must be > 0.
* The so called CIGAR string in the sixth field contains important mapping information. If these are unavailable, 
indicated by a "*", we regard the read as unmapped.


Read pairs are only considered if both sequences are mapped. Side note: The example
above complies with none of these criteria (indicated in red), hence this read is mapped.
For further processing, the CIGAR string, the mapping positions, the nucleotide sequence
(10th field), and the qualities in ASCII representation (11th field) are extracted.
The CIGAR string includes different operations to describe the mapped read. Regions
with alignment matches are denoted by "M", which does not necessarily mean that the
nucleotides agree with the reference sequence. An alignment match includes sequence
matches and mismatches. Deletions are denoted by the character "D", insertions by "I",
and soft clipped regions, which can not be mapped, are indicated by the character "S". The
preceding numbers in front of the respective characters denote the amount of consecutive
nucleotides to which this operation applies. In the example above, the CIGAR string
contains 94M6S, i.e. 94 positions are matching the alignment and six positions could not
be mapped. For construction of the ultimate read, only matching regions are considered.

## Merging
The matching sequence pair is merged to one single read for the count routine. However,
these reads may contain regions or single positions which are ignored, because they do
not fulfil the following requirements.
First of all, a symbol in the sequence other than the four bases A, C, G, and T is invalid.
This is also the case for nucleotides where the quality score is smaller than the given
threshold.
The pair of sequences can either be distinctly mapped to the reference such that there is a
region in between, which is not covered by the reads.
In this case, the nucleotides of the matching positions of both aligned sequences are simply
added to the merged read. In case indels occur in between the matching regions, they
are also treated as not covered regions.

However, the mapped sequences may also overlap. For the overlapping region, different cases are contemplated to 
decide for the incorporated nucleotide. If the two nucleotides are equal, which may be either the wild type of the
reference sequence or a mutation, the respective nucleotide is added to the read. If they
differ, but one of the sequences contains the wild type, the other one is considered as sequencing
error and the wild type is chosen for the read. This applies also, if the one wild
type nucleotide is given and the second symbol is invalid. If the two nucleotides disagree,
but both differ from the reference (including being invalid) the position is neglected.

## Counting
The constructed new read is evaluated according to the given dimension parameter.
In the case of single site variation (dimension = 1), the respective nucleotide count is
incremented by one for each incorporated position. After iterating through all reads, the
nucleotide occurrences (#) of all positions are written into the given output file, seperated
by tabulator:
```
pos #A #C #G #T
```
For the 2d case, co-occurring nucleotides for each position pair are counted. The resulting
output contains the frequencies in the following format:
```
pos1 pos2 #AA #AC #AG #AT #CA #CC ... #GG #GT #TA #TC #TG #TT
```
As example, the column "#AG" contains the number of reads which span the two respective
positions and comprise the nucleotide "A" at the first position and the nucleotide "G"
at the second position.
Similarly, the reads for a dimension of 3 are counted for each triplet of positions and nucleotide
co-occurrences and written into the output file:
```
pos1 pos2 pos3 #AAA #AAC #AAG #AAT #ACA ... #TGT #TTA #TTC #TTG #TTT
```
