# ped_from_relatedness
Attempt to construct small PEDs from relatedness statstics

## INSTALL

This package requires python3. To install via pip run:

    pip3 install git+git://github.com/gantzgraf/ped_from_relatedness.git

Alternatively, download or clone this repository and run:
    
    python3 setup.py

## INSTRUCTIONS

    usage: ped_from_relatedness RELATEDNESS VCF [options]

Two input files are required:
    
    RELATEDNESS - output from vcftools --relatedness or --relatedness2 algorithm
    VCF         - a bgzip compressed and indexed VCF file containing X 
                  chromosome variants for your cohort

Ideally your relatedness data should have been generated using variants that 
are relatively common (e.g. MAF >= 1%) in the general population. You may, 
however, modify the --family_cutoff, --first_degree_cutoff and --dup_cutoff
arguments to change the default relatedness cutoffs to suit your data. 

You may also wish to use the --min_gq and/or --pass options to this script to 
filter out poor quality X chromosome calls for gender checking. If your VCF is 
large, you may specify a maximum number of X chromosome calls to check using 
the --max_x_vars option.

## EXAMPLES

To attempt to generate a PED file indicating familial relationships from a set
of relatedness data using default cutoffs:

    ped_from_relatedness out.relatedness variants.vcf.gz > inferred.ped

As above, but data are from the vcftools --relatedness2 algorithm:

    ped_from_relatedness out.relatedness2 variants.vcf.gz -r2 > inferred.ped

As above, but only check the first 2000 non-PAR variants on the X chromosome 
when inferring gender:

    ped_from_relatedness out.relatedness2 variants.vcf.gz -r2 -n 2000 > inferred.ped

Check an existing PED file for potential errors:

    ped_from_relatedness out.relatedness2 variants.vcf.gz -r2 -n 2000 -p test.ped > ped_check.tsv

For detailed options see below.

## USAGE

    usage: ped_from_relatedness RELATEDNESS VCF [options]

    Construct or check small pedigrees from relatedness data.

    positional arguments:
      RELATEDNESS           Input relatedness/relatedness2 file.
      VCF                   Input VCF file containing at a minimum, genotypes on
                            the X chromosome. Must be bgzip compressed and tabix
                            indexed.

    optional arguments:
      -h, --help            show this help message and exit
      -p PED, --ped PED     PED file. If provided, relationships in this PED file
                            will be checked and flagged if problematic rather than
                            outputting the default output.
      -r2, --relatedness2   Input is output of --relatedness2 command from
                            vcftools. Default is to assume input is from the
                            --relatedness option of vcftools. This adjusts the
                            default relatedness cutoffs accordingly.
      -f FAMILY_CUTOFF, --family_cutoff FAMILY_CUTOFF
                            Custom relatedness value cutoff to consider two
                            samples related. Depending on the type of data used
                            (e.g. WES vs WGS) you may need to tune this threshold
                            to find a sensible value. Default value is 0.125 or
                            0.08 if --relatedness2 flag is set.
      -1 FIRST_DEGREE_CUTOFF, --first_degree_cutoff FIRST_DEGREE_CUTOFF
                            Custom relatedness value cutoff to consider two
                            samples first degree relatives. Default value is 0.35
                            or 0.2 if --relatedness2 flag is set.
      -d DUP_CUTOFF, --dup_cutoff DUP_CUTOFF
                            Custom relatedness value cutoff to consider two
                            samples as potential duplicates. Default value is 0.8
                            or 0.45 if --relatedness2 flag is set.
      --duplicates_file DUPLICATES_FILE
                            Write any detected duplicates to this file.
      -x XX_CUTOFF, --xx_cutoff XX_CUTOFF
                            Custom cutoff for ratio het variants to total variants
                            for assigning a sample as XX rather than XY. By
                            default, KMeans is used to determine genders from the
                            ratios of X chromosome heterozygous to homozygous
                            variants in your samples if you have at least 50
                            samples in your VCF, otherwise an arbitrary cutoff of
                            0.3 is used. This option overrides the use of KMeans
                            and the default arbitrary cutoff
      -a ASSEMBLY, --assembly ASSEMBLY
                            Assembly VCF input. Used for defining PAR region when
                            determining genders. Default=hg38.
      -m MIN_GQ, --min_gq MIN_GQ
                            Minimum genotype quality. Genotype calls with a GQ
                            below this value will be ignored. Default=20.
      -n MAX_X_VARS, --max_x_vars MAX_X_VARS
                            Maximum number of variants from the X chromosome to
                            test when checking genders. Default behaviour is to
                            check all variants in the non-pseudoautosomal region.
                            Use this option if your VCF is large and checking X
                            chromosome variants takes a longer than desired.
      --pass_filters        Only use variants with 'PASS' in the filter field for
                            inferring gender.

## AUTHOR

Written by David A. Parry at the University of Edinburgh. 

