TraCS: Transplant pair Comparison System
VERSION: Version 1.0 (May 30, 2013)

GENERAL USAGE NOTES
---------------------
PURPOSE: Directly compare two individual DNA samples for single nucleotide polymorphism (SNP) differences with an alloreactivity potential.

REQUIRED INPUT PARAMETERS:
	infile:    Multisample vcf file
        outdir:    Directory into which all of the data will reside
        direction: Direction of a potential immunological reaction (HVG, GVH, DI)
        db:        Path to the human database provided by Annovar

OPTIONAL INPUT PARAMETERS:
	ref:      Tab-delimited text file containing sample names (two names per line) to be compared
        bed:      Analyze a set of sites on the basis of a BED file.  Only the first three columns (chrom, chromStart and chromEnd) 		  	  are required. The BED file should have a header line.
        min:      Integer representing the minimum read coverage depth (default: 10)
        max:      Integer representing the maximum read coverage depth (default: 500)
        ver:      Annovar genome database version (default: hg18)
        threads:  The number of threads to be used for multiprocessing
        syn:      Flag to print out synonymous mutations in the GFF3 file. (default: OFF)
        keep:     Flag to keep all of the output files generated. (default: OFF)

OUTPUT FILES:
        counts file:     Tab-delimited text file containing the sample pair and count for each functional genotype.
        gff3 file:       Tab-delimited text file.
--------------------------------------------------------------------------------------------------------------------------------------


INSTALLATION
-----------------

TraCS requires Perl v5.8.8 or higher, UNIX command line interface, VCFtools v0.1.9, and ANNOVAR v2012-03-08.
Please include in your bash profile the paths to Perl and VCFtools. For example:

PATH=$PATH:/path/to/vcftools_0.1.9


CONTACT INFORMATION
----------------------
Contact:  Nihar Sheth
Phone:  (804) 827-0951
Email:  nsheth@vcu.edu
Website:  https://github.com/nisheth/TraCS

