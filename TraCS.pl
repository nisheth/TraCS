#!/usr/local/bin/perl
# TraCS: Transplant pair Comparison System
#    VERSION: Version 1 (March 11, 2013)
#    PURPOSE: Directly compare two individual DNA samples for single nucleotide
#             polymorphism (SNP) differences with an alloreactivity potential.
#
#    INPUT PARAMETERS:
#             infile:    Multisample vcf file
#             outdir:    Directory into which all of the data will reside
#             direction: Direction of a potential immunological reaction
#             db:        Path to the human database provided by Annovar
#
#    OPTIONAL INPUT PARAMETERS:
#             ref:      Tab-delimited text file containing sample
#                       names (two names per line) to be compared
#             bed:      Analyze a set of sites on the basis of a BED file.
#                       Only the first three columns (chrom, chromStart and chromEnd)
#                       are required. The BED file should have a header line.
#             min:      Integer representing the minimum read coverage depth
#             max:      Integer representing the maximum read coverage depth
#             ver:      Annovar database version
#             threads:  The number of threads to be used for multiprocessing
#             syn:      Flag to print out synonymous mutations in the GFF3 file.
#             keep:     Flag to keep all of the output files generated.
#
#    OUTPUT FILES:
#             .counts file:     Tab-delimited text file containing the sample pair
#                               and count for each functional genotype.
#             .gff3 file:       Tab-delimited text file.

############## LIBRARIES AND PRAGMAS ################

  use strict;
  use warnings;
  use File::Basename;
  use Getopt::Long;
  use threads;
  use Thread::Semaphore;

#################### VARIABLES ######################

my $infile = "";             # Multisample vcf file
my $outdir = "";             # Directory to store output files
my $db = "";                 # Path to ANNOVAR human database
my $ver = "hg18";            # Annovar reference genome build version
my $ref = "";                # File containing a list of sample names for reference
my $bed = "";                # Optional bed file
my $direction = "";          # Optional direction of change
my $help;                    # Readme help file
my $threads = 1;             # Number of threads to use
my $minDP = 10;              # Minimum cutoff for read coverage depth
my $maxDP = 500;             # Maximum cutoff for read coverage depth
my $syn = "";                # Switch to print synonymous mutations
my $keep = "";               # Switch to keep all output files
my $outfile = "";            # Output filename
my $pair;                    # Pair name
my $cmd = "";                # Statement to be executed on command line
my $pass = 0;                # Switch to determine if program needs to delete file(s)
my $bedfile;                 # Basename of bed file
my @location;                # Contains descriptive name of chromosome location(s)
my %ID;                      # Contains rsID numbers

################### MAIN PROGRAM ####################
#    Filter multisample vcf file and separate individual samples.
#    Transform genotypes into Annovar accepted input.
#    Perform a pairwise comparison of SNPs for each individual sample.
#    Annotate the SNPs using Annovar and sort based on function.

my $result = GetOptions('f|file=s' => \$infile,
                        'o|outdir=s' => \$outdir,
                        'a|annovar=s' => \$db,
                        'r|ref=s' => \$ref,
                        'min=i' => \$minDP,
                        'max=i' => \$maxDP,
                        'b|bed=s' => \$bed,
                        'd|direction=s' => \$direction,
                        'v|ver=s' => \$ver,
			't|threads=i' => \$threads,
                        'syn' => \$syn,
			'keep' => \$keep,
                        'h|help|man' => \$help
);

exit help_message() unless ($infile =~ m/\.vcf$/ && -d $outdir && -d $db && ($direction eq "GVH" || $direction eq "HVG" || $direction eq "DI"));

my $indir = $outdir;
$indir =~ s/\///;

if($bed ne ""){
        print "BED file provided for region-specific sample comparisons: $bed\n";
        $bedfile = basename($bed);
        @location = split(/\./, $bedfile);
        $outfile = basename($infile);
        $outfile =~ s/\.vcf/\.$location[0]/;
        $outfile = $outdir.$outfile;
        $cmd = "vcftools --vcf $infile --bed $bed --recode --out $outfile";     ## filter vcf file for regions listed in bed file
        runSystemCommand($cmd);
        $infile = $outfile.".recode.vcf";
}

print "Filtering $infile\n";
my $common_pos;
$infile = filterVCF($infile, $minDP, $maxDP, $bed);

my ($pref, $iref) = getSamplePairs($infile, $ref);
my @pairs = @$pref;
my @indiv = @$iref;

print "Separating individuals from $infile\n";
threadProcess('separateIndividual', $infile, $bed, \@indiv);

$indir = $outdir;
$indir =~ s/\/$//;

print "Transforming genotype data\n";
my @infiles = <$indir/*.MDR.SNP.filtered.recode.vcf>;
threadProcess('genotypeFix', \@infiles);
if(!$keep){ unlink @infiles; }

print "Comparing samples\n";
@infiles = <$indir/*.MDR.SNP.filtered.genotyped.recode.vcf>;
if($ref ne ""){ compareSamples(\@infiles, $ref); }
else{ compareSamples(\@infiles); }

print "Removing data which is the same as the reference\n";
threadProcess('removeRef', \@infiles);

print "Annotating SNPs\n";
@infiles = <$indir/*.MDR.SNP.filtered.genotyped.noRef.recode.vcf>;
threadProcess('SNPannotation', $db, $ver, \@infiles);
if(!$keep){ unlink @infiles; }

print "Adding annotations to comparison files\n";
my @diff = <$indir/*.diff.sites_in_files>;
@infiles = <$indir/*.exonic_variant_function>;
annotateDiff(\@infiles, \@diff, $direction);

my @avinput = <$indir/*.avinput*>;
my @idx = <$indir/*.vcfidx>;
my @log = <$indir/*.log>;
my @diff_indv = <$indir/*.diff.indv_in_files>;
my @SNP = <$indir/*.SNP.*>;
if(!$keep){
    unlink @avinput;
    unlink @diff;
    unlink @infiles;
    unlink @idx;
    unlink @log;
    unlink @diff_indv;
    unlink @SNP;
}

#################### SUBROUTINES ####################

##### FILTER MULTISAMPLE VCF FILE ####
#     Filter the multisample vcf based on the minimum
#     and maximum read coverage depth and remove
#     insertions/deletions.  Determine the number of
#     SNP positions common among all samples.

 sub filterVCF{
    my ($infile, $minDP, $maxDP, $bed) = (@_);
    my ($outfile, $geno_out, $common);

    $outfile = basename($infile);
    $outfile = $outdir.$outfile;
    if($bed ne ""){ $outfile =~ s/\.recode\.vcf/\.SNP\.filtered/; }
    else{ $outfile =~ s/\.vcf/\.SNP\.filtered/; }

    $cmd = "vcftools --vcf $infile --minDP $minDP --maxDP $maxDP --remove-indels --recode --out $outfile";
    runSystemCommand($cmd);                                                     ## filter for only SNPs and min/max coverage depth

    return $outfile.".recode.vcf";
 }

##### GET SAMPLE PAIRS ####
#     Obtain samples pairs to be analyzed from
#     either a reference file or multisample file.

 sub getSamplePairs{
    my($infile, $ref) = (@_);
    my($line, @cols, @pairs, @indiv, $outer, $inner);


    if($ref ne ""){
        open RF, "$ref" or die "Error opening $ref for reading\n";
        while($line = <RF>){
                chomp $line;
                @cols = split("\t", $line);
                if(scalar @cols != 2){ die "ERROR: Reference file should be tab-delimited with two columns.\n\tThe first column should contain a sample name to be compared to the sample name in the second column.\n"; }
		push @pairs, $cols[0]."-".$cols[1];
		if(!grep(/\A$cols[0]/, @indiv)){ push @indiv, $cols[0]; }
		if(!grep(/\A$cols[1]/, @indiv)){ push @indiv, $cols[1]; }
        }
        print "Sample pairs to be analyzed:  ", join("\t", @pairs), "\n";
    }
    else{
        print "WARNING: No reference file provided. Pairwise comparison of all samples.\n";
        @indiv = getIndividualNames($infile);                        # Get a list of the sample names from multisample vcf file
        for($outer = 1; $outer < (scalar @indiv); $outer++){
                for($inner = 0; $inner < $outer; $inner++){
                        push(@pairs, $indiv[$inner]."-".$indiv[$outer]);
                }
        }
        print "Sample pairs to be analyzed:  ", join("\t",@pairs), "\n";
    }
    return (\@pairs, \@indiv);
 }

#### GET INDIVIUAL NAMES ####
#    Read in the multisample .vcf file.
#    Extract the individual sample names.

 sub getIndividualNames{
   my $infile = shift;
   open FH, "$infile" or die "Error opening $infile for reading\n";

   my ($iline, @icols, @individual);
   while($iline = <FH>){
        if($iline =~ /CHROM/g){
                chomp $iline;
                @icols = split(/\t/, $iline);
                my $size = @icols;
                for(my $i = 9; $i < $size; $i++){
                        $individual[$i - 9] = $icols[$i];
                }
                close FH;
                return @individual;
        }
   }
 }

##### SEPARATE INDIVIDUALS ####
#     Separates individuals from the multisample
#     vcf file and removes any missing data denoted
#     by ./.

 sub separateIndividual{
    my ($infile, $bed, $sem, $sample) = (@_);
    my ($outfile, $new_infile, $out);
    
    $outfile = $infile;
    if($bed ne ""){ $outfile =~ s/.+\.$location[0]/$outdir$sample\.$location[0]/; $outfile =~ s/\.recode\.vcf//; }
    else{ $outfile =~ s/.+\.vcf/$outdir$sample\.SNP\.filtered/; }
    $out = $outfile.".recode.vcf";
    if(!-e $out){
	$cmd = "vcftools --vcf $infile --indv $sample --recode --out $outfile";
	runSystemCommand($cmd);                                         ## separate samples into individual files

        $new_infile = $outfile.".recode.vcf";
        $outfile = $new_infile;
        $outfile =~ s/\.SNP/\.MDR\.SNP/;
        $outfile =~ s/\.recode\.vcf//;
        $cmd = "vcftools --vcf $new_infile --geno 1 --recode --out $outfile";       ## filter out missing data from sample file
        runSystemCommand($cmd);
    }
    $sem->up; # Used for threading
 }
 
#### THREAD PROCESS ####
#    Create a thread for a subroutine

 sub threadProcess{
   my @params = (@_);
   my ($sub, $data_array, $sem, @threads);
   $sub = shift @params;
   $data_array = pop @params;
   $sem = Thread::Semaphore->new($threads); # max allowable threads
   @threads = map {
        $sem->down;
        if(scalar @params > 0){ threads->create($sub, @params, $sem, $_) }
        else{ threads->create($sub, $sem, $_) }
   } @$data_array;
   $_->join for @threads; # Finish remaining threads in use
 }

##### GENOTYPE FIX ####
#     Transform the genotype and alternate allele to
#     to an Annovar-accepted format

 sub genotypeFix{
    my ($sem, $infile) = (@_);
    my $output = $infile; $output =~ s/.filtered./.filtered.genotyped./;
    open FH, "$infile" or die "Error opening $infile for reading\n";
    open OFH, ">$output" or die "Error opening $output for writing\n";

    my ($line, @cols, @geno, @num, @allele);
    while($line = <FH>){
        chomp $line;
        if($line =~ /^chr/){
                @cols = split(/\t/, $line);
                @geno = split(/:/, $cols[9]);
                if($cols[2] ne "."){ $ID{$cols[0]}{$cols[1]} = $cols[3]; }      ## Store the rsID number
                if($cols[4] =~ /,/g){                                           ## Alt allele column has more than one allele denoted
                        if($geno[0] eq "0/0"){ $cols[4] = $cols[3]; }           ## If genotype is 0/0, set alt allele same as ref allele
                        @num = split(/\//, $geno[0]);
                        if(($geno[0] =~ /0/g) || ($num[0] == $num[1])){         ## If the genotype info contains 0 or is homozygous
                                if($geno[0] =~/1/g){                            ## If genotype is 0/1, 1/0, or 1/1
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[0];                  ## set the 2nd allele as alt allele
                                }
                                elsif($geno[0] =~/2/g){                         ## If genotype is 0/2, 2/0, or 2/2
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[1];                  ## set the 2nd allele as alt allele
                                        $geno[0] =~ s/2/1/g;                    ## replace the 2 in the genotype with a 1
                                        $cols[9] =~ s/^[0-9]\/[0-9]/$geno[0]/g; ## replace the original genotype with new genotype
                                }
                                elsif($geno[0] =~/3/g){                         ## if genotype is 0/3, 3/0, or 3/3
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[2];                  ## set the 3rd allele as alt allele
                                        $geno[0] =~ s/3/1/g;                    ## replace the 3 in the genotype with a 1
                                        $cols[9] =~ s/^[0-9]\/[0-9]/$geno[0]/g; ## replace the original genotype with new genotype
                                }
                        }
                        else{                                                   ## If the genotype info is not homozygous or heterozygous to reference
                                if($geno[0] eq "1/3" || $geno[0] eq "3/1"){     ## if genotype is 1/3, 3/1
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[0].",".$allele[2];   ## set the 1st and 3rd allele as alt allele
                                        $geno[0] =~ s/3/0/g;                    ## replace the 3 in the genotype with a 0
                                        $cols[9] =~ s/^[0-9]\/[0-9]/$geno[0]/g; ## replace the original genotype with new genotype
                                }
                                elsif($geno[0] eq "1/2" || $geno[0] eq "2/1"){  ## if genotype is 1/2, 2/1
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[0].",".$allele[1];   ## set the 2nd and 3rd allele as alt allele
                                        $geno[0] =~ s/2/0/g;                    ## replace the 2 in the genotype with a 0
                                        $cols[9] =~ s/^[0-9]\/[0-9]/$geno[0]/g; ## replace the original genotype with new genotype
                                }
                                elsif($geno[0] eq "2/3" || $geno[0] eq "3/2"){  ## if genotype is 2/3, 3/2
                                        @allele = split(/,/, $cols[4]);
                                        $cols[4] = $allele[1].",".$allele[2];   ## set the 2nd and 3rd allele as alt allele
                                        $geno[0] =~ s/2/0/g;                    ## replace the 2 in the genotype with a 0
                                        $geno[0] =~ s/3/1/g;                    ## replace the 3 in the genotype with a 1
                                        $cols[9] =~ s/^[0-9]\/[0-9]/$geno[0]/g; ## replace the original genotype with new genotype
                                }
                        }
                }
                else{                                                     ## Alt allele column has only 1 allele
                        if($geno[0] eq "0/0"){ $cols[4] = $cols[3]; }           ## If genotype is 0/0, set alt allele same as ref allele
                }
                print OFH "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$cols[7]\t$cols[8]\t$cols[9]\n";
        }
        else{ print OFH "$line\n"; }
    }
    close FH;
    close OFH;

    $sem->up; # Used for threading
 }

#### COMPARE SAMPLES ####
#    If no reference file is provided, all samples will
#       undergo a pairwise comparison.
#    Otherwise, use the list of sample pair names provided.
#    Perform a pairwise comparison of the samples listed
#       in the reference file.

 sub compareSamples{
   my $files = shift;
   my $list = shift;
   my @files = @$files;
   my ($sline, @scols, @file, $base, @oArray, @rArray, @pArray, $outer, $inner);

   if(defined($list)){
        open RF, "$list" or die "Error opening $list for reading\n";
        while($sline = <RF>){
                chomp $sline;
                @scols = split("\t", $sline);
                if(!grep(/\/$scols[0]/, @files)){ die "Error: \"".$scols[0]."\" is not a sample in the multisample vcf file.\n"; }
                elsif(!grep(/\/$scols[1]/, @files)){ die "Error: \"".$scols[1]."\" is not a sample in the multisample vcf file.\n"; }
                if(grep(/\/$scols[0]/, @files)){ @file = grep(/\/$scols[0]/, @files); push @rArray, $file[0]; }
                if(grep(/\/$scols[1]/, @files)){ @file = grep(/\/$scols[1]/, @files); push @oArray, $file[0]; }
                push @pArray, ($scols[0]."-".$scols[1]);
        }
        close RF;
   }
   ## Compare each reference file with each non-reference file.
   if(!defined($list)){
        @pArray = ();
        for($outer = 1; $outer < (scalar @files); $outer++){
                for($inner = 0; $inner < $outer; $inner++){
                        if(!grep/$files[$inner]/, @oArray){ push @oArray, $files[$inner]; }
                }
                threadProcess('vcfDiff', $files[$outer], \@pArray, \@oArray);
        }
   }
   else{ foreach my $ref(@rArray){ threadProcess('vcfDiff', $ref, \@pArray, \@oArray); } }
 }
 
#### VCF DIFF ####
#    Based on the number of paramters:
#       create an output file for the comparison.
#       run VCFtools to compare 2 samples.

 sub vcfDiff{
    my ($ref, $pairList, $sem, $other) = (@_);
    my (@pairs, $output, $sample, @sub, $cmd);
    @pairs = @$pairList;
    if(scalar @pairs != 0){
        foreach(@pairs){
                @sub = split("-", $_);
                if($ref =~ /\/$sub[0]\./ && $other =~ /\/$sub[1]\./){
                        $output = $other;
                        $output =~ s/$sub[1]/$_/;
                        $cmd = "vcftools --vcf $ref --diff $other --out $output";
                        runSystemCommand($cmd);
                }
        }
    }
    elsif($ref =~ /\/(.+)\.SNP/){
        $output = $other;
        $sample = $1;
        @sub = split(/\./, $sample);
        $output =~ s/(?<=\/)/$sub[0]-/;
        $cmd = "vcftools --vcf $ref --diff $other --out $output";
        runSystemCommand($cmd);
    }
    $sem->up; # Used for threading
 }

#### REMOVE REFERENCE ####
#    Remove positions where the genotype is the
#       as the reference denoted as "0/0".

 sub removeRef{
   my ($sem, $infile) = (@_);
   my $output = $infile; $output =~ s/.genotyped./.genotyped.noRef./;
   open FH, "$infile" or die "Error opening $infile for reading\n";
   open OFH, ">$output" or die "Error opening $output for writing\n";

   my($line, @cols, @geno);
   while($line = <FH>){
        chomp $line;
        if($line =~ /^chr/){
                @cols = split(/\t/, $line);
                if($cols[9] !~ /^0\/0/){ print OFH "$line\n"; }
        }
        else{ print OFH "$line\n"; }
   }
   close FH;
   close OFH;
   $sem->up; # Used for threading
 }

#### SNP ANNOTATION ####
#    Convert sample.SNP.filtered.genotyped.recode.vcf
#       files to an Annovar-accepted format.
#    Annotate the files using Annovar and the provided
#       human database.

 sub SNPannotation{
   my($humandb, $ver, $sem, $infile) = (@_);

   ## Convert files into an .avinput file
   my $avinput = $infile.".avinput";
   $cmd = "convert2annovar.pl -format vcf4 -allallele $infile > $avinput";
   runSystemCommand($cmd);

   ## Annotate the .avinput file
   $cmd = "annotate_variation.pl -buildver $ver $avinput $humandb >& /dev/null";
   runSystemCommand($cmd);
   
   $sem->up; # Used for threading
 }

#### ANNOTATE DIFFERENCE ####
#    Store the contents of the exonic_variant files in a hash.
#    Add annotation data from the exonic_variant hash to the
#    matching chromosomes and positions in the difference file.

 sub annotateDiff{
    my($exfiles, $diff, $change) = (@_);
    my @exfiles = @$exfiles;
    my @diff = @$diff;
    my %annoHash = processExonicFile(@exfiles);
    my ($output, $pos_output);
    
    if(defined($bedfile)){
        $output = $outdir."Samples.$location[0].$change.counts.txt";
        open OFH, ">$output" or die "Error opening $output for writing\n";
        print OFH "Pair\tLocation\tSynonymous\tNonsynonymous\tCons\tNonCons\tStop\tUnk\n";

        $pos_output = $outdir."Samples.$location[0].$change.positions.txt";
    }
    else{
        $output = $outdir."Samples.$change.counts.txt";
        open OFH, ">$output" or die "Error opening $output for writing\n";
        print OFH "Pair\tSynonymous\tNonsynonymous\tCons\tNonCons\tStop\tUnk\n";

        $pos_output = $outdir."Samples.$change.positions.txt";
    }

    threadProcess('addAnnoDir', \%annoHash, $direction, \@diff );
 }

#### PROCESS EXONIC FILE ####
#    Store annotation information based on sample,
#    chromosome, start position, and zygosity.

 sub processExonicFile{
        my @infile = @_;
        my ($sub, @sub, $line, @cols, %Fxn);
        foreach my $infile(@infile){
                open FH, "$infile" or die "Error opening $infile for reading\n";
                $sub = basename($infile);
                @sub = split(/\./, $sub);
                while($line = <FH>){
                        chomp $line;
                        @cols = split(/\t/, $line);
                        $Fxn{$sub[0]."_".$cols[3]}{$cols[4]}{$cols[7]} = "$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[6]\t$cols[7]\t$cols[8]";
                }
                close FH;
        }
        return %Fxn;
 }

#### ADD ANNOTATION AND DIRECTION ####
#    For each SNP chromosome and position, add
#    Annovar-annotation based on the direction
#    of the change and count the number of
#    SNPs based on the type of mutation.

 sub addAnnoDir{
    my ($annoHash_ref, $direction, $sem, $infile) = (@_);
    my($s1, $s2, $line, @cols, @anno1, @anno2, $s1allele, $s2allele, $ref, @s1allele, %counts);
    my %annoHash = %$annoHash_ref;

    open FH, "$infile" or die "Error opening $infile for reading\n";

    $pair = getPairName($infile);
    if($pair =~ /(.+)-(.+)/){ $s1 = $1; $s2 = $2; }

    <FH>;
    while($line = <FH>){
        chomp $line;
        @cols = split(/\t/, $line);
        if($cols[2] eq "B"){       ## if data is present for both samples
           if($direction eq "GVH"){ %counts = directionalChange($s1.'_'.$cols[0], $s2.'_'.$cols[0], $cols[1], $cols[3], $cols[4], $cols[5], \%counts, $annoHash_ref); }
           elsif($direction eq "HVG"){ %counts = directionalChange($s2.'_'.$cols[0], $s1.'_'.$cols[0], $cols[1], $cols[3], $cols[5], $cols[4], \%counts, $annoHash_ref); }
           elsif($direction eq "DI"){ %counts = directionalChange($s1.'_'.$cols[0], $s2.'_'.$cols[0], $cols[1], $cols[3], $cols[4], $cols[5], \%counts, $annoHash_ref); }
        }
    }
    my @Fxn = ('Syn', 'Nonsyn', 'Cons', 'NonCons', 'Stop', 'Unk');
    my($chr, $fxn, $total);
    if($bed ne ""){
        print OFH "$pair\t$location[0]\t";
        foreach $chr(sort keys %counts){
                foreach $fxn(@Fxn){
                        if(!exists($counts{$chr}{$fxn})){ $counts{$chr}{$fxn} = 0; }
                        print OFH "$counts{$chr}{$fxn}\t";
                }
                print OFH "\n";
        }
    }
    else{
          print OFH "$pair\t";
          foreach $fxn(@Fxn){
                $total = 0;
                foreach $chr(sort keys %counts){
                        if(!exists($counts{$chr}{$fxn})){ $counts{$chr}{$fxn} = 0; }
                        $total = $total + $counts{$chr}{$fxn};
                }
                print OFH "$total\t";
          }
          print OFH "\n";
    }
    $sem->up; # Used for threading
 }

#### DIRECTIONAL CHANGE ####
#    Given sample 1 and sample 2 alleles, if the sample 1
#    allele(s) has the potential to produce a protein that
#    is foreign in comparison with sample 2, then there is
#    a directional change.

 sub directionalChange{
    my ($s1_chr, $s2_chr, $pos, $ref, $s1allele, $s2allele, $chk, $anno) = (@_);
    my %chk = %$chk;    my %anno = %$anno;
    my (@anno1, @anno2, @s1allele);

    if($s1allele !~ /,/){                       ## if s1 has 1 allele
        if(exists $anno{$s1_chr}{$pos}{$s1allele}){
                @anno1 = split(/\t/, $anno{$s1_chr}{$pos}{$s1allele});
                if($anno1[6] eq 'hom'){ $s1allele = $s1allele.','.$s1allele; }
                elsif($anno1[6] eq 'het'){ $s1allele = $s1allele.','.$ref; }
        }
        else{ $s1allele = $ref.','.$ref; }      ## s1 does not exist in hash therefore it is the same as the reference
    }
    else{ $s1allele = $s1allele; }              ## if s1 has 2 alleles
    if($s2allele !~ /,/){                       ## if s2 has 1 allele
        if(exists $anno{$s2_chr}{$pos}{$s2allele}){
                @anno2 = split(/\t/, $anno{$s2_chr}{$pos}{$s2allele});
                if($anno2[6] eq 'hom'){ $s2allele = $s2allele.','.$s2allele; }
                elsif($anno2[6] eq 'het'){ $s2allele = $s2allele.','.$ref; }
        }
        else{ $s2allele = $ref.','.$ref; }      ## s2 does not exist in hash therefore it is the same as the reference
    }
    else{ $s2allele = $s2allele; }              ## if s2 has 2 alleles
    $ref = $ref.','.$ref;                       ## make reference allele have 2 alleles
    @s1allele = split(/,/, $s1allele);

    if($s1allele ne $s2allele){                 ## s1 and s2 are SNPs to each other
        if($direction eq "DI"){
                if($s1allele eq $ref){ %chk = chkAAchangetoRef($s2_chr, $pos, $s2allele, $s1allele, \%chk, $anno, 0); }
                elsif($s2allele eq $ref){ %chk = chkAAchangetoRef($s1_chr, $pos, $s1allele, $s2allele, \%chk, $anno, 1); }
                else{ %chk = chkAAchange($s1_chr, $s2_chr, $pos, $s1allele, $s2allele, \%chk, $anno); }
        }
        elsif($s2allele !~ /$s1allele[0]/ || $s2allele !~ /$s1allele[1]/){
                if($s1allele eq $ref){ %chk = chkAAchangetoRef($s2_chr, $pos, $s2allele, $s1allele, \%chk, $anno, 0); }
                elsif($s2allele eq $ref){ %chk = chkAAchangetoRef($s1_chr, $pos, $s1allele, $s2allele, \%chk, $anno, 1); }
                else{ %chk = chkAAchange($s1_chr, $s2_chr, $pos, $s1allele, $s2allele, \%chk, $anno); }
        }
    }
    return %chk;
 }

#### CHECK AMINO ACID CHANGE TO REFERENCE ####
#    If the one of the sample alleles is the same
#    as the reference and the other is different,
#    then check the type of mutation that occurs

 sub chkAAchangetoRef{
    my ($s1_chr, $pos, $s1allele, $s2allele, $hash, $anno, $flip) = (@_);
    my %hash = %$hash;  my %anno = %$anno;
    my (@s1_1_anno, @s1_2_anno, $s2_aa, $s2_1_aa, $s1_aa, $s1_1_aa, $s1_2_aa, $s1_1_allele, $s1_2_allele);
    my %aa = (  'A' => 'non-polar', 'V' => 'non-polar', 'L' => 'non-polar',
                'I' => 'non-polar', 'M' => 'non-polar', 'G' => 'non-polar',
                'F' => 'aromatic',  'W' => 'aromatic',  'P' => 'aromatic',
                'Y' => 'aromatic',  'T' => 'polar',     'C' => 'polar',
                'N' => 'polar',     'Q' => 'polar',     'S' => 'polar',
                'K' => 'basic',     'R' => 'basic',     'H' => 'basic',
                'D' => 'acidic',    'E' => 'acidic' );

    ($s1_1_allele, $s1_2_allele) = split(/,/, $s1allele);
    if($anno{$s1_chr}{$pos}{$s1_1_allele}){
        @s1_1_anno = split(/\t/, $anno{$s1_chr}{$pos}{$s1_1_allele});

        if($s1_1_anno[0] eq "unknown"){           ## if the position is annotated as an unknown mutation, count as Unk
                $hash{$s1_1_anno[2]}{'Unk'}++;
                if($flip == 1){ printPos('unknown', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                else{ printPos('unknown', $s1_chr, $pos, $s2allele, $s1allele, $s2_aa, $s1_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
        }
        else{
                if($anno{$s1_chr}{$pos}{$s1_1_allele}){
                        if($s1_1_anno[1] =~ /p\.([A-Z])\d+([A-Z])/){ $s2_1_aa = $1; $s1_1_aa = $2; $s1_2_aa = $1; }
                        $s1_aa = $s1_1_aa.$s1_2_aa;
                }
                if($anno{$s1_chr}{$pos}{$s1_2_allele}){
                        @s1_2_anno = split(/\t/, $anno{$s1_chr}{$pos}{$s1_2_allele});
                        if($s1_2_anno[1] =~ /p\.[A-Z]\d+([A-Z])/){ $s1_2_aa = $1; }
                        $s1_aa = $s1_1_aa.$s1_2_aa;
                }
                ## create two allele amino acids
                $s2_aa = $s2_1_aa.$s2_1_aa;

                if($s2_1_aa eq 'X' || $s1_1_aa eq 'X' || $s1_2_aa eq 'X'){        ## if either s1 allele 1 aa, s1 allele 2 aa, or s2 aa is an X, count as Stop
                        $hash{$s1_1_anno[2]}{'Stop'}++;
                        if($flip == 1){ printPos('stop', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                        else{ printPos('stop', $s1_chr, $pos, $s2allele, $s1allele, $s2_aa, $s1_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                }
                elsif(($s2_1_aa eq $s1_1_aa) && ($s2_1_aa eq $s1_2_aa)){    ## if s2 aa, allele 1 aa, and allele 2 aa are all equal, then count as Syn
                        $hash{$s1_1_anno[2]}{'Syn'}++;
                        if($syn ne ""){
                                if($flip == 1){ printPos('synonymous', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                                else{ printPos('synonymous', $s1_chr, $pos, $s2allele, $s1allele, $s2_aa, $s1_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                        }
                }
                else{       ## if s1 alleles are not Unk, Syn, or Stop, count as Nonsyn
                        $hash{$s1_1_anno[2]}{'Nonsyn'}++;
                        if($aa{$s2_1_aa} eq $aa{$s1_1_aa} && $aa{$s2_1_aa} eq $aa{$s1_2_aa}){       ## if both 1st and 2nd alleles are conservative changes, count as Cons
                                $hash{$s1_1_anno[2]}{'Cons'}++;
                                if($flip == 1){ printPos('conservative', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                                else{ printPos('conservative', $s1_chr, $pos, $s2allele, $s1allele, $s2_aa, $s1_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                        }
                        else{   ## count as NonCons
                                $hash{$s1_1_anno[2]}{'NonCons'}++;
                                if($flip == 1){ printPos('nonconservative', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                                else{ printPos('nonconservative', $s1_chr, $pos, $s2allele, $s1allele, $s2_aa, $s1_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                        }
                }
        }
    }
    return %hash;
 }

#### CHECK AMINO ACID CHANGE ####
#    Check the type of mutation that occurs
#    at each amino acid for each allele

 sub chkAAchange{
    my ($s1_chr, $s2_chr, $pos, $s1allele, $s2allele, $hash, $anno) = (@_);
    my %hash = %$hash; my %anno = %$anno;
    my (@s1_1_anno, @s1_2_anno, @s2_1_anno, @s2_2_anno, $s1_1_allele, $s1_2_allele, $s2_1_allele, $s2_2_allele);
    my ($s1_aa, $s1_1_aa, $s1_2_aa, $s2_aa, $s2_1_aa, $s2_2_aa);
    my %aa = (  'A' => 'non-polar', 'V' => 'non-polar', 'L' => 'non-polar',
                'I' => 'non-polar', 'M' => 'non-polar', 'G' => 'non-polar',
                'F' => 'aromatic',  'W' => 'aromatic',  'P' => 'aromatic',
                'Y' => 'aromatic',  'T' => 'polar',     'C' => 'polar',
                'N' => 'polar',     'Q' => 'polar',     'S' => 'polar',
                'K' => 'basic',     'R' => 'basic',     'H' => 'basic',
                'D' => 'acidic',    'E' => 'acidic' );

    ($s1_1_allele, $s1_2_allele) = split(/,/, $s1allele);
    ($s2_1_allele, $s2_2_allele) = split(/,/, $s2allele);
    if(exists $anno{$s1_chr}{$pos}{$s1_1_allele} && exists $anno{$s2_chr}{$pos}{$s2_1_allele}){
        @s1_1_anno = split(/\t/, $anno{$s1_chr}{$pos}{$s1_1_allele});
        @s2_1_anno = split(/\t/, $anno{$s2_chr}{$pos}{$s2_1_allele});

        if($s1_1_anno[0] eq 'unknown' || $s2_1_anno[0] eq 'unknown'){   ## if either s1 or s2 allele has an unknown mutation, count as Unk
                $hash{$s1_1_anno[2]}{'Unk'}++;
                printPos('unknown', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele});
        }
        else{
                if($s1_1_anno[1] =~ /p\.([A-Z])\d+([A-Z])/){ $s1_1_aa = $2; $s1_2_aa = $1; }
                if(exists $anno{$s1_chr}{$pos}{$s1_2_allele}){
                        @s1_2_anno = split(/\t/,$anno{$s1_chr}{$pos}{$s1_2_allele});
                        if($s1_2_anno[1] =~ /p\.[A-Z]\d+([A-Z])/){ $s1_2_aa = $1; }
                }
                $s1_aa = $s1_1_aa.$s1_2_aa;     ## create two allele amino acids

                if($s2_1_anno[1] =~ /p\.([A-Z])\d+([A-Z])/){ $s2_1_aa = $2; $s2_2_aa = $1; }
                if(exists $anno{$s2_chr}{$pos}{$s2_2_allele}){
                        @s2_2_anno = split(/\t/,$anno{$s2_chr}{$pos}{$s2_2_allele});
                        if($s2_2_anno[1] =~ /p\.[A-Z]\d+([A-Z])/){ $s2_2_aa = $1; }
                }
                $s2_aa = $s2_1_aa.$s2_2_aa;     ## create two allele amino acids

                if($s1_1_aa eq 'X' || $s1_2_aa eq 'X' || $s2_1_aa eq 'X' || $s2_2_aa eq 'X'){
                        $hash{$s1_1_anno[2]}{'Stop'}++;
                        printPos('stop', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele});
                }
                elsif($s1_1_aa eq $s2_1_aa && $s1_1_aa eq $s2_2_aa && $s1_2_aa eq $s2_1_aa && $s1_2_aa eq $s2_2_aa){
                        $hash{$s1_1_anno[2]}{'Syn'}++;
                        if($syn ne ""){ printPos('synonymous', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele}); }
                }
                else{
                        $hash{$s1_1_anno[2]}{'Nonsyn'}++;
                        ## if all alleles are conservative changes, count as Cons
                        if($aa{$s2_1_aa} eq $aa{$s1_1_aa} && $aa{$s2_1_aa} eq $aa{$s1_2_aa} && $aa{$s2_2_aa} eq $aa{$s1_1_aa} && $aa{$s2_2_aa} eq $aa{$s1_2_aa}){
                                $hash{$s1_1_anno[2]}{'Cons'}++;
                                printPos('conservative', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele});
                        }
                        else{   ## count as NonCons
                                $hash{$s1_1_anno[2]}{'NonCons'}++;
                                printPos('nonconservative', $s1_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno{$s1_chr}{$pos}{$s1_1_allele});
                        }
                }
        }
    }
    return %hash;
 }

#### PRINT POSITIONS ####
#    Print the chromosome, mutation type, position,
#    and a list of attributes for each SNP to a
#    separate file for further analysis.

 sub printPos{
    my($type, $s_chr, $pos, $s1allele, $s2allele, $s1_aa, $s2_aa, $anno) = (@_);
    my($filename, $sample, $chr, $color, @anno, @anno_1, $v, $aa, $pos_0);
    ($sample, $chr) = split(/_/, $s_chr);

    if($bed ne ""){ $filename = $outdir.$pair.".$location[0]".".$direction".".gff3"; }
    else{ $filename = $outdir.$pair.".$direction".".gff3"; }

    if(-e $filename && $pass == 0){ unlink $filename or warn "Could not unlink $filename"; $pass = 1; }     ## remove an existing file if 1st time running
    if(-e $filename){ open PFH, ">>$filename" or die "Error opening $filename for appending\n"; }
    else{
        open PFH, ">$filename" or die "Error opening $filename for writing\n";
        if($bed ne ""){ print PFH "##gff-version 3\n#track name=".$pair.".$location[0]".".$direction"." itemRgb=On\n"; }
        else{ print PFH "##gff-version 3\n#track name=".$pair.".$direction"." itemRgb=On\n"; }
    }

    if($type eq 'conservative'){ $color = "#E8950A"; }          ## gold
    elsif($type eq 'nonconservative'){ $color = "#0000FF"; }    ## blue
    elsif($type eq 'stop'){ $color = "#FF0000"; }               ## red
    elsif($type eq 'unknown'){ $color = "#00FF00"; }            ## green
    elsif($type eq 'synonymous'){ $color = "#6D008C"; }         ## purple

    if($s1allele =~ /,/){ $s1allele =~ s/,//; }
    if($s2allele =~ /,/){ $s2allele =~ s/,//; }
    $v = $s2allele.'->'.$s1allele;

    $pos_0 = $pos - 1;

    if($type eq 'unknown'){
        if(exists $ID{$chr}{$pos}){
                print PFH "$chr\t.\t$type\t$pos_0\t$pos\t.\t.\t.\tVar=$v;AA=unknown;gene=unknown;ref=$ID{$chr}{$pos};color=$color\n";
        }
        else{ print PFH "$chr\t.\t$type\t$pos_0\t$pos\t.\t.\t.\tVar=$v;AA=unknown;gene=unknown;ref=.;color=$color\n"; }
    }
    else{
        $aa = $s2_aa.'->'.$s1_aa;
        @anno = split(/\t/, $anno);
        @anno_1 = split(/:/, $anno[1]);
        if(exists $ID{$chr}{$pos}){
                print PFH "$chr\t.\t$type\t$pos_0\t$pos\t.\t.\t.\tVar=$v;AA=$aa;gene=$anno_1[0];ref=$ID{$chr}{$pos};color=$color\n";
        }
        else{ print PFH "$chr\t.\t$type\t$pos_0\t$pos\t.\t.\t.\tVar=$v;AA=$aa;gene=$anno_1[0];ref=.;color=$color\n"; }
    }
    close PFH;
 }

#### GET PAIR NAME ####
#    Extract the name of the pair
#    being analyzed.

 sub getPairName{
    my $infile = shift;
    $infile = basename($infile);
    my @pair = split(/\./, $infile);
    return $pair[0];
 }

#### HELP MESSAGE ####
#    Prints a help message to the screen

  sub help_message{
     print "
        Usage:
                Masterpipeline.pl -f <filename> -o <output Directory> -a <path/to/annovar/humandb/> [option] [parameter]

                Required:
                        -f, --file		multisample vcf filename
                	-o, --outdir		directory for output files
                	-a, --annovar		path/to/ANNOVAR/humandb/
                	-d, --direction	        direction of the change (GVH, HVG, or DI)

                Optional:
                	-r, --ref               Tab-delimited file containing a list of sample pair names
                                                to use for comparison (ex: sample1      sample2).
                	-min			minimum cut off for read coverage depth (default: 10)
                	-max			maximum cut off for read coverage depth (default: 500)
                	-b, --bed		bed filename
                	-t, --threads		number of threads to use
                        -v, --ver               Annovar database version (default: hg18)
                        --syn                   Switch to print synonymous mutations in gff3 output file(s)
                                                (default: off)
                        --keep                  Switch to keep all files generated by the program
                                                (default: off)
                	-h, --help, --man	brief help message
                \n";
  }

#### RUN SYSTEM COMMAND ####
#    Executes command line commands
#    Prints the command to the screen

  sub runSystemCommand{
     my $cmd = shift;
#     print "CMD - $cmd\n";
#     print "\n";
     `$cmd`;
     return;
  }
