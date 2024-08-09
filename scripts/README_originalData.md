
* six_cell_lines_minimal(A549).csv

downloaded via emailed link from Oleksandra Fanari <fanari.o@northeastern.edu>
on 20240808: 

The filtering here is very minimal as we discussed today: pvalue<0.001 and a mmIVT<10  (5 of U-to-C mismatch) to exclude SNV. 
You can further filter this data to include only those positions with a mm.Direct >10.

If you want to narrow down to only high-confidence sites then you may increase the mm.Direct to 30 or 40%.


Columns meaning:
Annotation: gene name
chr: chromosome
position: index of modified U in hg38 coordinates
kmer: 5-mer with the modified U in the center
T_Direct: number of T in the direct sample for the considered position
C_Direct: same as before but for C
T_IVT: number of T in the IVT library
C_IVT: same for C
N_reads_Direct: number of reads for the position (A,C,T,G, ins and dels)
N_reads_IVT: same as before in the IVT library
expected.mm: you can disregard this, it's used to calculate the mm in the Direct/IVT library
mm.Direct: % of U-to-C mismatch detected by nanopore for that position. The higher the better as more reads are found to carry psi instead of a u.
mm.IVT:% of U-to-C mismatch in the IVT library. The lower the better as it means the unmodified RNA probably has an SNV.
p.value.Direct: statistical significance of that position being called a psi modification.
mm.DirectMINUSmm.IVT: mm.Direct-mm.IVT percentage of mismatch
total_Direct: number of T + number of C for that position in the Direct library
total_IVT: same as before for IVT
psi: number of psi sites found on that transcript
ACP: AnnotationChrPosition
TRUE: internal computations, you can disregard
which.IVT: the IVT data can come from the pairedIVT or from a panIVT composed by all the 6 cell lines to supplement for low number of reads in the in vitro trascribed library
CP:ChrPosition

