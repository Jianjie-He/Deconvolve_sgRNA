# Deconvolve_sgRNA
This is a Python script to deconvolve superimposed sequencing chromatograms of sgRNAs. It works well in chromatograms with one or two overlapped sgRNAs but loses effectiveness with three overlapped sgRNAs and malfunctions with four or more sgRNAs. For chromatograms with three overlapped sgRNAs, we recommend decomposing them manually with the Cartesian product. Alternatively, to deal with this problem caused by the inconsistent peak position of each base, you can set the base cut-off ratio below 0.10. This script is handy in sequencing chromatograms of sgRNAs generated in pooled CRISPR library transformation.
