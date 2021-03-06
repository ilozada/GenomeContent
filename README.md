# __GenomeContent__


A program that estimates global statistics and sequence-based estimators of genome features: exons, introns and intergenic regions

A brief introduction to the program
-----------------

`GenomeContent` was written in Perl to calculate global statistics and sequence-based estimators of genome features in six major steps, as shown in Figure 1. First, the processing of gene annotations focuses on identifying coordinates from protein-coding gene (CDS), while the filtering process focuses on checking the quality of intron annotations. As described in next sections, the coordinates derived from both proceses are taken as the “reference gene sets” for introns, exons and intergenic regions to directly estimate several statistic descriptors, such as size, density and number. Then, the “reference gene sets” are projected onto the genome sequence in both strands, so that the nucleotide contents of each genome-feature are calculated according to the definitions described in a section below. Finally, all statistic descriptors obtained with the program are provided as text files, fasta formats and exploratory figures (as shown in Figure 1). `GenomeContent` runs on an entire genome in few minutes or hours, depending on genome size and the number of annotated genes.

---

### Irma Lozada-Chávez

*Please send your questions, comments and bug reports to: : irma.lozada.chavez@gmail.com*

If you use `GenomeContent` for your research, please cite the following paper:

- Lozada-Chávez I, Stadler,P.F., and Prohaska, S.J. (2018) [__Genome-wide features of introns are evolutionary decoupled among themselves and from genome size throughout Eukarya__](https://www.biorxiv.org/content/early/2018/03/18/283549). __under peer-review__. doi: https://doi.org/10.1101/283549

___

Table of Contents
-----------------

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Program options](#program-options)
4. [Results: short overview](#results-short-overview)
5. [Manual](#Manual)
6. [Trouble shooting](#trouble-shooting)
6. [Description](#description)


## Installation

`GenomeContent` v1.0 is written in Perl language. So, it can be easily modified by the user. Additional documentation along the code is provided to facilitate these modifications. The program is released as a compressed archive file (`GenomeContent.tar.gz`) for download. Extracting the files to your current directory will create the directory `GenomeContent`, containing the required files.

`Unzip`:
```bash
tar -zxvf GenomeContent.tar.gz
```

Perl version (and above) and the following libraries are required to run the program (linux only):
[Perl](https://www.perl.org/get.html), [Getopt::Long](http://search.cpan.org/~jv/Getopt-Long-2.49.1/lib/Getopt/Long.pm), [POSIX](http://search.cpan.org/~bingos/perl/ext/POSIX/lib/POSIX.pod), [Math::Complex](http://search.cpan.org/~zefram/Math-Complex-1.59/lib/Math/Complex.pm), [Chart::Gnuplot](http://search.cpan.org/dist/Chart-Gnuplot/lib/Chart/Gnuplot.pm), [GD::Graph](http://search.cpan.org/dist/GDGraph/Graph.pm).

## Quick start

**Input files** \
Two files are required by `GenomeContent`, the **genome sequence in fasta** format (with the file extension: **faa**, **fas**, **fna** or **fasta**), and the **annotation file in General Feature or Transfer Formats** (with the file extension: **gtf**, **gff** or **gff3**) where the coordinates of the protein-coding genes are described. There are two versions of the GFF file format in general use: GFF v2 (created by the Sanger Institute) and GFF v3 (created by the Sequence Ontology Project). Both versions have a number of differences to consider. [See complete description of `input` files](#input-files-details).

The `GenomeContent` program can be run in two different modes: **single mode** and **non-single mode**. The **single mode** is designed to perform on one single genome project, this is, one genome sequence with its corresponding annotation file. The **non-single mode** is designed to perform over several genome projects as one job, this is, two or more different annotation files located in one directory with their corresponding genome sequences. See complete description of the `program options` files [here](#program-options) and `output` files [here](#output-files).

`GenomeContent`: General and mandatory options:

  ```terminal
$ perl GenomeContent.pl

  MANDATORY OPTIONS
  -s|single <yes|no>          Mode to be run: single = yes or no. Single mode is for one genome project, while the non-single mode is for several genomes.
  -a|annotations <dir/file>   Path to and name of the annotation file in GFF or GTF format
  -g|genomes <dir/file>       Path to and name of the genome file in FASTA format
  -o|outputs <dir/file>       Path to and name of the directory to write the output files
  -d|database <name>          Name of the database where the annotation file was downloaded
  -s|status <genome|assembly> Status of the sequenced genome project: genome or assembly
  -e|species <name>           The name of the species as written in the input files


  OPTIONAL PARAMETERS
  -me <integer>               Minimum nucleotide length for an exon. Default: 15
  -mi <integer>               Minimum nucleotide length for an intron. Default: 15
  -h|help                     Displays help
```


- **SINGLE MODE** \
  *Syntaxis:*

```terminal
perl GenomeContent.pl -s yes -a <annotations> -g <genomes> -o <outfiles> -d <database> -t <status> -e <species> [optional parameters]
```

- **NON-SINGLE MODE** \
  *Syntaxis:*

```terminal
  perl GenomeContent.pl  -s no -a <annotations> -g <genomes> -o <outfiles> [optional parameters]
```

**Test `GenomeContent`: quick tutorial** \
Examples `infiles` and `outfiles` to test `GenomeContent` program are located in the directory named [Examples](directory/where/infiles_test/are_stored/EXaMPLEs). Use these examples to run the program (after [installing](#installation) all prerequisite perl libraries) in:

`single mode` (one genome):

```terminal
perl GenomeContent.pl -s no -a <annotations> -g <genomes> -o <outfiles> [optional parameters]
```

`non-single mode` (several genomes):

```terminal
perl GenomeContent.pl -s y -a /home/user/volvox_carteri.genes.gff -g /home/user/volvox_carteri.assembly.fasta -o /home/user/outputs
```

#### Input files: details

Two files are required by `GenomeContent`, the **genome sequence in fasta** format (with the file extension: **faa**, **fas**, **fna** or **fasta**), and the **annotation file in General Feature or Transfer Formats** (with the file extension: **gtf**, **gff** or **gff3**) where the coordinates of the protein-coding genes are described. There are two versions of the GFF file format in general use: GFF v2 (created by the Sanger Institute) and GFF v3 (created by the Sequence Ontology Project). Both versions have a number of differences to consider.

**GFF2** represents one nesting level of features and only two-level feature hierarchies, such as: transcipt → exon, so that `GenomeContent` identifies in the “group field of hierarchies” the name and ID numbers of the “CDS” feature coordinates that come from the same gene. The GTF format is a refinement of GFF2 and is sometimes referred to as GFF2.5. PROBLEM IN ALTERNATIVE TRANSCRIPTS.

**GFF3** supports arbitrarily many hierarchical levels and  three-level feature hierarchies, such as: gene → transcript → exon. The top level is a feature of type "gene" which bundles up the gene's transcripts and regulatory elements. Beneath this level are one or more transcripts of type "mRNA". This level can also accommodate promoters and other cis-regulatory elements. At the third level are the components of the mRNA transcripts, most commonly CDS coding segments and UTRs.
This example shows how to represent a gene named "EDEN" which has three alternatively-spliced mRNA transcripts:

For example, it  and gives specific meanings to certain tags in the attributes field.

**NOTE:** In principle,

Once the user has both files, it is important to make sure that these files fulfill three requirements:

**1.** Make sure that the name or ID of the chromosomes/scaffolds/contigs in the fasta format correspond to the same names located in the first column of the gene annotation file. *Examples*:
Fasta file 1:

```terminal
>17 dna:chromosome chromosome:GRCh37:17:1:81195210:1 REF
AAGCTTCTCACCCTGTTCCTGCATAGATAATTGCATGACAATTGCCTTGTCCCTGCTGAA
TGTGCTCTGGGGTCTCTGGGGTCTCACCCACGACCAACTCCCTGGGCCTGGCACCAGGGA
```

Gene annotation file 1:

```terminal
17  protein_coding  CDS     17991338        17991401    .    +       0       gene_id "ENSG00000108591"; transcript_id "ENST00000580929"; exon_number "1"; gene_name "DRG2"; gene_biotype "protein_coding"; transcript_name "DRG2-011"; protein_id "ENSP00000462060";
```

Fasta file 2:

```terminal
>scaffold_1
CAGAACGGGAAACAGAAGAAAATCGTGTGAAGACGAAATAATATCGCTGGCCAGTACGGCCGCGGTACATTCAACATGTA
AAACTTTCTCATGTCTGTCTATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCC
```

Annotation file 2:

```terminal
scaffold_1		phytozome8_0	CDS	8410627	8410668	.	-	0	ID=PAC:23125614.CDS.1; Parent=PAC:23125614;pacid=23125614
```

**NOTE:** In principle, both files should correspond to the same released version of the genome project. In case that the user would like to compare different versions of these files, just be aware that some protein-coding gene coordinates could not be mapped correctly into the genome sequence if changes in the contigs, scaffolds or chromosomes occur along the versions.

**2.** If the size of the genome sequence is bigger than 2.5 Gigabytes and the user is using a single standard laptop (2.6GHz, 8GB RAM), the user might consider to divide the genome sequence in smaller files. For example,the pinus genome ≈ 22.0 Gbs will need ~60 Gbs to read and process the sequence. Thus, the user is advised to divide the genome sequence by chromosomes or parts of similar size.

**3.** Make sure that the names of the sequence and the annotation files have the following structure:

```terminal
species.database.status.otheruserinfo.fasta  and  species.database.status.otheruserinfo.gff

<species>     Name of the species, e.g.: homo, homo_sapiens, hsa
<database>    Name of the database where the annotation file was downloaded, e.g.: ensembl, phytozome, jgi, vectorbase, denovo, etc.
              The program makes a particular distinction between the “ensemble” database and the rest. Make sure that you provide the correct database name.

<status>      One of the only two status recognized for a genome project:
              genome   = if complete until chromosome assignments
              assembly = if assembled at the conting and/or scaffold level
              The program requires any of these two options to be correctly provided for the single mode.
              The option should also match the corresponding information provided in the name of the annotation file.
<fasta>       Extension of the fasta file, e.g.: .fasta, .faa, .fas, .fna
<annotation>  Extension of the fasta file, e.g.: .gff, .gtf
```

Example for an assembled genome:

```bash
volvox_carteri.phytozome.assembly.gff  and  volvox_carteri.phytozome.assembly.fasta
```

Example for a complete genome in a single file or several files:

```bash
homo_sapiens.ensembl.genome.gff  and  homo_sapiens.ensembl.genome.fasta         or
                                      homo_sapiens.ensembl.chromosome.1.fasta
                                      homo_sapiens.ensembl.chromosome.2.fasta
									  ...
                                      homo.sapiens.ensembl.chromosome.22.fasta
                                      homo.sapiens.ensembl.chromosome.X.fasta
                                      homo.sapiens.ensembl.chromosome.Y.fasta
```

**NOTE:** `GenomeContent` requires that all names of the files keep the above structure, regardless the mode you choose to run the pipeline, this is for a single or several genome projects. The structure of the name is particularly important when `GenomeContent` performs over several genome projects in one call.


#### Output files

In the specified folder, `-o`, several output-files are generated. A description of each output file is now provided with examples for [the green algae Volvox carteri genome project v2.0 (Prochnik et al., 2010)](https://www.ncbi.nlm.nih.gov/pubmed/20616280) from the [Phytozome database version 10.0](https://phytozome.jgi.doe.gov/pz/portal.html). It is important to note that the figures are just provided as a guide to explore the data. Instead, the information provided in the text files is, thus, more suitable to be used to plot figures in a higher quality with other plotting programs.

`GenomeContent` creates several `output` files. See the complete description of each `output` file [here](#description-output-files). Complete list of the output files:

- *outfile.coordinates.txt*
- *outfile.size.txt*
- *outfile.[exons|introns|intergenics].fasta*
- *outfile.filtered.exons.fasta | outfile.filtered.introns.fasta*
- *outfile.[exons|introns|intergenics].content.txt*
- *outfile.exonsintronsintergenics.content.txt*
- *outfile.all.sizes.[exons|introns|intergenics].distribution.txt*
- *outfile.all.sizes.[exons|introns|intergenics].distribution.ps*
- *outfile.quartiles.sizes.[exons|introns|intergenics].distribution.txt*
- *outfile.quartiles.sizes.[exons|introns|intergenics].distribution.ps*
- *outfile.order.[exons|introns|intergenics].distribution.txt*
- *outfile.exons-introns.order.sizes.ps*
- *outfile.quartiles.order.[exons|introns].distribution.txt*
- *outfile.exons-introns.order.sizes.ps*
- *outfile.quartiles.order.[exons|introns].distribution.txt*
- *outfile.[quartiles.]order.introns-exon-intron.distribution.txt*
- *outfile.[quartiles.]exons-introns.definition.sizes.ps*
- *outfile.density.distribution.[exons|introns].genome.ps*
- *total.outfile.[exons|introns|intergenics]_ranges.table.txt*
- *outfile.ranges.sizes.distribution.genome.ps*
- *total.outfile.intron_types.table.txt | total.outfile.intron_types.strands.table.txt*
- *total.outfile.[exons|introns|intergenics]_ranges_atgc.table.txt*
- *total.outfile.summary.genome.txt*
- *total.outfile.[exons|introns|intergenics].table.txt*
- *total.outfile.[exons|introns|intergenics].strands.table.txt*
- *total.outfile.[exons|introns|intergenics].chromosomes.table.txt*


## Program options

The `GenomeContent` program can be run in two different modes: **single mode** and **non-single mode**. The **single mode** is designed to perform on one single genome project, this is, one genome sequence with its corresponding annotation file. The **non-single mode** is designed to perform over several genome projects as one job, this is, two or more different annotation files located in one directory with their corresponding genome sequences.

General options:

```terminal
MANDATORY OPTIONS
-s|single <yes|no>          Mode to be run: single = yes or no. Single mode is for one genome project, while the non-single mode is for several genomes.
-a|annotations <dir/file>   Path to and name of the annotation file in GFF or GTF format
-g|genomes <dir/file>       Path to and name of the genome file in FASTA format
-o|outputs <dir/file>       Path to and name of the directory to write the output files
-d|database <name>          Name of the database where the annotation file was downloaded
-s|status <genome|assembly> Status of the sequenced genome project: genome or assembly
-e|species <name>           The name of the species as written in the input files


OPTIONAL PARAMETERS
-me <integer>               Minimum nucleotide length for an exon. Default: 15
-mi <integer>               Minimum nucleotide length for an intron. Default: 15
-h|help                     Displays help
```

- **SINGLE MODE**

```terminal
perl GenomeContent  -s yes -a <annotations>  -g <genomes>  -o <outfiles>
                    -d <database>  -t <status>  -e <species> [optional parameters]


-a  <annotations>   Path to and name of the annotation file in GFF or GTF format.
                    [REQUIRED] e.g., -a /home/user/volvox_carteri.genes.gff

-g  <genomes>       Path to and name of the genome file in FASTA format.
                    [REQUIRED] e.g., -g /home/user/volvox_carteri.assembly.fasta

-o  <outfiles>      Path to and name of the directory to write the output files.
                    [REQUIRED] e.g., -o /home/user/outputs

-d  <database>      Name of the database where the annotation file was downloaded.
                    [REQUIRED] e.g., -d phytozome

-t  <status>        Status of the sequenced genome.
                    [REQUIRED] e.g., -t assembly

-e  <species>       The name of the species as written in the input files.
                    [REQUIRED] e.g., -e volvox_carteri

-me <minexon>       The minimum nucleotide size required for exons. Default = 15 nts.
                    [OPTIONAL] e.g., -me 50

-mi <minintron>     The minimum nucleotide size required for introns. Default = 15 nts.
                    [OPTIONAL] e.g., -mi 25
```


*Syntaxis:*

```termial
perl GenomeContent.pl  -s yes  -a /home/user/volvox_carteri.genes.gff -g /home/user/volvox_carteri.assembly.fasta -o /home/user/outputs -d phytozome -t assembly -e volvox_carteri
```

- **NON-SINGLE MODE**

```terminal
perl GenomeContent.pl  -s no -a <annotations> -g <genomes> -o <outfiles> [optional parameters]


-a  <annotations>   Path to the annotation files in GFF or GTF format.
                    [REQUIRED] e.g., -a /home/user/annotations

-g  <genomes>       Path to and name of the genome file in FASTA format.
                    [REQUIRED] e.g., -g /home/user/genomes

-o  <outfiles>      Path to and name of the directory to write the output files.
                    [REQUIRED] e.g., -o /home/user/outputs

-me <minexon>       The minimum nucleotide size required for exons. Default = 15 nts.
                    [OPTIONAL] e.g., -me 50

-mi <minintron>     The minimum nucleotide size required for introns. Default = 15 nts.
                    [OPTIONAL] e.g., -mi 25
```

*Syntaxis:*

```terminal
perl GenomeContent.pl  -s y  -a /home/user/volvox_carteri.genes.gff -g /home/user/volvox_carteri.assembly.fasta  -o /home/user/outputs
```

## Results: short overview

In the specified folder, -o, several output-files are generated. A description of each output file is now provided with examples for the green algae [Volvox carteri genome project v2.0](https://www.ncbi.nlm.nih.gov/pubmed/20616280) (Prochnik et al., 2010, PMID: 20616280) from the [Phytozome database version 10.0](https://phytozome.jgi.doe.gov/pz/portal.html). It is important to note that the figures are just provided as a guide to explore the data. Instead, the information provided in the text files is, thus, more suitable to be used to plot figures in a higher quality with other plotting programs.

#### Description output files

- **outfile.coordinates.txt** \
It reports the coordinates and lengths for introns and exons of the protein coding genes as well as for the intergenic regions within each sequence. The quartile where the length of every single feature falls is also provided; quartiles are estimated from the total feature population (e.g., the total number of introns or exons or intergenics).

*Example of the outfile.coordinates.txt:*

```terminal
SEQ_NAME     ID_CDS               FEATURE       START       END     STRAND SIZE    QUARTILE
scaffold_1   vcn:PACid:23124778   CDS_exon_1    14008323    14008400   +   78      Q0:<93
scaffold_1   vcn:PACid:23124778   CDS_exon_2    14008491    14008610   +   120     Q1-Q2:93-144
scaffold_1   vcn:PACid:23124778   CDS_exon_3    14008699    14008813   +   115     Q1-Q2:93-144
scaffold_1   vcn:PACid:23124778   CDS_exon_4    14008980    14009401   +   422     Q4:>231
scaffold_1   vcn:PACid:23124778   CDS_intron_1  14008401    14008490   +   90      Q0:<205
scaffold_1   vcn:PACid:23124778   CDS_intron_2  14008611    14008698   +   88      Q0:<205
scaffold_1   vcn:PACid:23124778   CDS_intron_3  14008814    14008979   +   166     Q0:<205
scaffold_1   vcn               	  intergenic_1  1           31757      +   31757   Q4:>16491
scaffold_1   vcn               	  intergenic_2  32997       41019      +   8023    Q2-Q3:6848-16491
scaffold_1   vcn                  intergenic_3  47801       54377      +   6577    Q1-Q2:2269-6848
scaffold_1   vcn               	  intergenic_4  61967       62657      +   691     Q0:<2269
```

- **outfile.size.txt** \
It contains the nucleotide size of each fasta sequence included in the file of the genome sequence.

*Example of the outfile.size.txt:*

```terminal
scaffold_1      14152940
scaffold_100    83517
scaffold_1005   2386
scaffold_1006   2384
```

- **outfile.[exons|introns|intergenics].fasta** \
Each file provides the nucleotide sequence of all exons, introns and intergenic regions in fasta format, respectively. The description of the fasta format in the three files contain the following information:

```terminal
>gene_ID | location= seq_name : start-end : strand : feature_number : intron_type(only for introns) | organism= name_species | cds= length(nts) : exon_length= length_nts (% from total CDS size)  | A= length_nts (% from total nts size) :  T= length_nts (% from total nts size) : G= length_nts (% from total nts size) : C= length_nts (% from total nts size) : AT= length_nts (% from total nts size) : N= length_nts (% from total nts size)
```

*Example of the outfile.exons.fasta:*

```terminal
>vcn:PACid:23124778|location=scaffold_1:14008323-14008400:+:CDS_exon_1 |organism=volvox_carteri|cds=1079(nts):exonic=735(nts)|exon_length=78(7.23%)|A=19(24.36%):T=17(21.79%):G=19(24.36%):C=23(29.49%):AT=36(46.15%):GC=42(53.85%):N=0(0.00%)
ATGCAACATTCCATTAGCAAAAGTGTGGAGCTCCGACGGCGCTGCCCAGGGCCTTGTTCCACTAGCTACACGATTACG
...
>vcn:PACid:23124778|location=scaffold_1:14008980-14009401:+:CDS_exon_4 |organism=volvox_carteri|cds=1079(nts):exonic=735(nts)|exon_length=422(39.11%)|A=96(22.75%): T=79(18.72%):G=132(31.28%):C=115(27.25%):AT=175(41.47%):GC=247(58.53%):N=0(0.00%)
ATATGCGACCCTCACCCAACGTGGAAATACAGCGTCAAGAGGAAAGCGCACCAGGCATGTACCGGTACTATGAGCGCATCCAAATCAACAGCAGTGGCAGCTACAGGGGTTACATCAGCCAGCCGGCAGCTGCTGCGACGTCGGTAGTCAAACCGGTGTCGCCTTCCCTCAACCCCGCGTTAGTGGGCGCTCTGGTTGTGCTGGGCGGCTACGTGACCCTCACTGCGGCATTCAACCGCAATTACGCTCTTACCAACTACGCGGAAAGCAAGCGCTGGCTGCTGATGGCGTTGTGGCCTCTGTTGTTTCTGTTCTCCCCAAAGTTTAGGGAGCAGTTTGCAGCAGCGATTCGCGGAGAGAGGGTGGGGCTGAAGCGCGAAGGGACAAGGGAGGGGGACAGCAACGCTGGCAGCAGCGAATAG
```

*Example of the outfile.introns.fasta:*

```terminal
>vcn:PACid:23124778|location=scaffold_1:14008401-14008490:+:CDS_intron_1:3n |organism=volvox_carteri|cds=1079(nts):intronic=344(nts)|intron_length=90(8.34%)|A=22(24.44%):T=26(28.89%):G=23(25.56%):C=19(21.11%):AT=48(53.33%):GC=42(46.67%):N=0(0.00%)
GTGGGGATGCCCTCCAAGGGACCTTTCAGGTCGGCCACTAGGCATGATTGCAAACAGTTTGTAATCCTTTGATGTCAAAATTATTTGCAG
...
>vcn:PACid:23124778|location=scaffold_1:14008814-14008979:+:CDS_intron_3:3n+1 |organism=volvox_carteri|cds=1079(nts):intronic=344(nts)|intron_length=166(15.38%)|A=24(14.46%):T=41(24.70%):G=44(26.51%):C=57(34.34%):AT=65(39.16%):GC=101(60.84%):N=0(0.00%)
GTGTGTCAGCGTAGTCGTAGACGCGGTGAAGGGCCGGAAAATCTGTCCGTCACTGCCCCGACCAGCCTTTGGTCGGTATGCACACCGTTCGCCGTTTTGTTCGGGCGGAATCCCCCCTAATTCCCCATCCTGCACCCTTTGCTGTCCGTTCCCGCTCCGCCTGAAG
```

*Example of the outfile.intergenics.fasta:*

```terminal
>vcn:scaffold_1:+:1-31757:intergenic_1(31757 nts) |organism=volvox_carteri|length=31757(0.22%) |A=6784(21.36%):T=6612(20.82%):G=8249(25.98%):C=7324(23.06%) :AT=13396(42.18%):GC=15573(49.04%):N=2788(8.78%)
CAGAACGGGAAACAGAAGAAAATCGTGTGAAGACGAAATAATATCGCTGGCCAGTACGGCCGCGGTACATTCAACATGTAAAACTTTCTCATGTCTGTCTATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCCATCTACCCATCTATCTATCTACCCATCTATTTAGCTATCTATCTAGCTATCTGTCTACCCATCTATCTATCTACCCATCTATTTAGCTATCTATCTAGCTATCTATCTAGCTATATAAAACCCTATTA...
```

- **outfile.filtered.exons.fasta | outfile.filtered.introns.fasta** \
Each file provides the nucleotide sequence of all filtered exons and introns in fasta format, respectively. The default minimum length of an exon or an intron is 15 nts. The options `-me and -mi` should be used to change this parameter (e.g., -me 50 and -mi 25 will filter exons and introns smaller than 50 nts and 25 nts, respectively). The description of the fasta format in both files contain the same information that  `outfile.filtered.[exons|introns].fasta`.

- **outfile.[exons|introns|intergenics].content.txt** \
Each file provides the coordinates of all overlapped intervals in the sequence (i.e., chromosome or scaffold or contig) harboring exons, introns and intergenic regions, respectively. The information is provided like fasta format with the following information:

```terminal
>seq_name | seq_length | feature: total_number : number_after_overlap | feature_content_nts : feature_content_Mbs : feature_content_fraction_from_total_sequence_size
```

*Example of the outfile.exons.content.txt:*

```terminal
>vcn:scaffold_1|14152940_nts|exons:13859:13611|3052731_nts:3.05_Mbs:21.57%
25078:25145 25228:25484 26359:26481 27986:28152 31758:31885 31992:32121 32265:32355 32860:32996 35205:35383 35596:35737 36273:36353 36650:36739 37091:37159 37496:37542
... 14109835:14110134 14110629:14110703 14115906:14116229 14116319:14116345
>vcn:scaffold_100|83517_nts|exons:4:4|1494_nts:0.00_Mbs:1.79%
3563:3809 3827:4323 18057:18356 18368:18817
>vcn:scaffold_1005|2386_nts|exons:1:1|957_nts:0.00_Mbs:40.11%
4:960
```

*Example of the outfile.introns.content.txt:*

```terminal
>vcn:scaffold_1|14152940_nts|introns:12020:11833|4595203_nts:4.60_Mbs:32.47%
25146:25227 25485:26358 26482:27985 31886:31991 32122:32264 32356:32859 35384:35595 35738:36272 36354:36649 36740:37090 37160:37495 37543:38022 38156:38643 38777:39731
... 4099459 14099786:14100242 14101339:14101523 14101760:14101956 14102172:14102254
>vcn:scaffold_100|83517_nts|introns:1:1|17_nts:0.00_Mbs:0.02%
3810:3826
>vcn:scaffold_1008|2376_nts|introns:10:10|523_nts:0.00_Mbs:22.01%
391:417 566:625 870:885 1034:1093 1190:1249 1294:1353 1450:1509 1606:1665 1762:1821 1918:1977
```

*Example of the outfile.intergenics.content.txt:*

```terminal
>vcn:scaffold_1|14152940_nts|intergenics:1767|6504910_nts:6.50_Mbs:45.96%
1:25077 28153:31757 32997:35204 40280:41019 47801:54377 61967:62657 64596:65905 72239:74169 77316:77489 83735:84129 ... 14110704:14115905 14116346:14152940
>vcn:scaffold_100|83517_nts|intergenics:3|81995_nts:0.08_Mbs:98.18%
1:3562 4324:18056 18818:83517
>vcn:scaffold_1005|2386_nts|intergenics:2|1429_nts:0.00_Mbs:59.89%
1:3 961:2386
>vcn:scaffold_1007|2383_nts|intergenics:1|2383_nts:0.00_Mbs:100%
1:2383
>vcn:scaffold_1008|2376_nts|intergenics:2|597_nts:0.00_Mbs:25.13%
1:293 2073:2376
```

- **outfile.features.content.distribution.genome.ps** \
The information represented by the image is two fold. On the left axis, it is plotted the fraction or proportion of nucleotides from the respective genomic feature (intron, exon or intergenic regions) that contributes to genome size (depicted in Mbs). On the right axis, the nucleotide composition (AT, GC, Ns) in percentage for each genomic feature is plotted from the total nucleotide size of the respective genomic feature (depicted in Mbs).


<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.features.content.distribution.genome.png">
</p>
[Figure Content]

- **outfile.all.sizes.[exons|introns|intergenics].distribution.txt** \
Each file provides the frequencies of occurrence of all lengths for exons or introns or intergenic sizes.

  *Example of the outfile.all.sizes.exons.distribution.genome.txt:*

```text
  EXON_SIZE FREQUENCY
  100 433
  101 469
  102 785
  103 503
  104 465
  105 817
```

  *Example of the outfile.all.sizes.introns.distribution.genome.txt:*

```text
  INTRON_SIZE FREQUENCY
  100 111
  101 104
  102 144
  103 122
  104 107
  105 130
```

  *Example of the outfile.all.sizes.exons.distribution.genome.txt:*

  ```text
  INTERGENIC_SIZE FREQUENCY
  100 1
  102 1
  103 1
  105 1
  106 2
  107 4
  ```

- **outfile.all.sizes.[exons|introns|intergenics].distribution.ps** \
Each image shows the frequencies of occurrence of all lengths for exons or introns or intergenic sizes.

    *Example of outfile.all.sizes.exons.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.all.sizes.exons.distribution.genome.png">
</p>
[Figure Content]

    *Example of outfile.all.sizes.introns.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.all.sizes.introns.distribution.genome.png">
</p>
[Figure Introns Distribution]

  *Example of outfile.all.sizes.intergenics.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.all.sizes.intergenics.distribution.genome.png">
</p>
[Figure Intergenics Distribution]

- **outfile.quartiles.sizes.[exons|introns|intergenics].distribution.txt** \
Each file provides the frequencies of those lengths for exons or introns or intergenic sizes that are located between the `lower fence`<sup>&dagger;</sup> and the `upper fence`<sup>&Dagger;</sup> of the population distribution. Thus, outliers (*i.e.*, *extreme values*) for each feature are filtered in these files.

    <sup>&dagger;: Lower fence = first quartile – (1.5 * (third quartile – first quartile)).</sup><br>
    <sup>&Dagger;: Upper fence = third quartile + (1.5 * (third quartile – first quartile)).</sup>

- **outfile.quartiles.sizes.[exons|introns|intergenics].distribution.ps** \
Each image shows the frequencies of those lengths for exons or introns or intergenic sizes that are located between the lower fence and the upper fence of the population distribution. Different measurements of mean size are also provided: the standard, the weighted, the quartile and the normalised average size. See section 5. “Table for the description of the headers in the output files" for the definitions of the different mean measurements.

  *Example of outfile.quartiles.sizes.exons.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.quartiles.sizes.exons.distribution.genome.png">
</p>
[Figure Exons Sizes]

*Example of outfile.quartiles.sizes.introns.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.quartiles.sizes.introns.distribution.genome.png">
</p>
[Figure Introns sizes]

*Example of outfile.quartiles.sizes.intergenics.distribution.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.quartiles.sizes.intergenics.distribution.genome.png">
</p>
[Figure Intergenics sizes]

- **outfile.order.[exons|introns|intergenics].distribution.txt** \
Each file provides some statistic estimators about the average size of exons and introns according to their position in the protein-coding gene. See section 5. “Table for the description of the headers in the output files" for the description of the variables.

*Example of the outfile.order.exons.distribution.txt:*

```terminal
EXON_POSITION   TOTAL_EXONS   AVERAGE_SIZE   FILTERED_SIZE   MIN_SIZE   MAX_SIZE   SIZE_Q_DOWN   	SIZE_Q1   		SIZE_Q2-MEDIAN   SIZE_Q3   SIZE_Q_UP   SIZE_LOWER_FENCE   SIZE_UPPER_FENCE   	SIZE_VARIANCE   SIZE_SD   SIZE_ERRORBAR_BEGIN   SIZE_ERRORBAR_END   SIZE_FILTER_VAR   	SIZE_FILTER_SD   SIZE_FILTER_EBAR_BEGIN   SIZE_FILTER_EBAR_END
exon_1   15285   281.86   207.26   15   7555   39   104   195   325   836   -227.5   656.5   
         118629.70  	344.43   -62.57   626.29   19587.36   139.95   67.31   347.21
exon_2   12826   235.90   156.65   15   6948   48   93   145   241   730   -129   463   111312.39
         333.64   -97.74   569.54   8666.17   93.09   63.56   249.74
exon_3   11494   213.70   144.43   15   9746   47   90   137   216   641   -99   405   86983.97
         294.93   -81.23   508.63   6358.22   79.74   64.69   224.17
```

*Example of the outfile.order.introns.distribution.txt:*

```terminal
INTRON_POSITION   TOTAL_INTRONS   AVERAGE_SIZE   FILTERED_SIZE   MIN_SIZE   MAX_SIZE   SIZE_Q_DOWN   	SIZE_Q1   SIZE_Q2-MEDIAN   SIZE_Q3   SIZE_Q_UP   SIZE_LOWER_FENCE   SIZE_UPPER_FENCE   	SIZE_VARIANCE   SIZE_SD   SIZE_ERRORBAR_BEGIN   SIZE_ERRORBAR_END   SIZE_FILTER_VAR  	SIZE_FILTER_SD   SIZE_FILTER_EBAR_BEGIN   SIZE_FILTER_EBAR_END
intron_1   12806   377.59   322.40   15   3886   65   144   295   487   983   -370.5   1001.5					118316.55   343.97   33.62   721.56   487   32.78   220.76   101.64   543.16
intron_2   11482   409.86   352.78   15   3589   73   187   331   521   1028   -314   1022   
			117066.73			342.15   67.71   752.01   479   79.59   219.04   133.74   571.82
intron_3   10270   417.50   360.91   15   3465   76   204   343   528   1029   -282   1014   
			112861.65			335.95   81.55   753.45   460   37.54   214.56   146.35   575.47
```

- **outfile.exons-introns.order.sizes.ps** \
The image shows the standard average size for all the exons and introns by their position in a protein-coding gene.

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.exons-introns.order.sizes.png">
</p>
[Figure Exon-Intron order sizes]

- **outfile.quartiles.order.[exons|introns].distribution.txt** \
Each file provides some statistic estimators from the feature population (exons or introns) located between the lower fence and the upper fence and according to their position in the protein-coding gene. See section 5. “Table for the description of the headers in the output files" for the description of the variables.

- **outfile.exons-introns.order.sizes.ps** \
The image shows the standard average size for those exons and introns located between the lower fence and the upper fence and according to their position in a protein-coding gene. The frequency of occurrences for each position category is also plotted.

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.quartiles.exons-introns.order.sizes.png">
</p>
[Figure Exon-Intron order sizes]

- **outfile.exons-introns.order.sizes.ps** \
The image shows the error bars of the average sizes for those intron and exon located between the lower fence and the upper fence and according to their position in a protein-coding gene.

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.errorbars.quartiles.exons-introns.order.sizes.png">
</p>
[Figure Exon-Intron order sizes]

- **outfile.[quartiles.]order.introns-exon-intron.distribution.txt** \
These files provide the lengths of the upstream and downstream introns surrounding each internal exon within a protein-coding gene. The file labeled with “quartiles” includes the data only for those exons and introns located between the lower fence and the upper fence of the size distribution.

*Example of the outfile.order.intron-exon-intron.distribution.txt:*

```terminal
SPECIE:ID				EXON_NUMBER	EXON_SIZE		INTRON-UP_SIZE	INTRON-DOWN_SIZE
vcn:PACid:23124778		2				120				90					88
vcn:PACid:23124778		3				115				88					166
vcn:PACid:23124780		2				1290			964					670
vcn:PACid:23124780		3				651				670					393
vcn:PACid:23124780		4				672				393					186
vcn:PACid:23124780		5				123				186					195
```

- **outfile.[quartiles.]exons-introns.definition.sizes.ps** \
These images display internal exons as a function of the nucleotide length of its upstream (x axis) and downstream (y axis) introns. Each point represents an unique exon. The user migth see one of two patterns: most exons are flanked by long introns (> 250 nts), or most exons are flanked by at least one short intron (≤ 250 nts). The 250 nts demarcate the experimentally determined 200- to 250-nt transition from cross-intron to cross-exon splice site recognition (Fox-Walsh et al., 2005, PMID: 16260721). The file name labeled with “quartiles” (down image) includes the data only for those exons and introns located between the lower fence and the upper fence.

*Example of outfile.exons-introns.definition.sizes.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.exons-introns.definition.sizes.png">
</p>
[Figure exons-introns definition sizes]

*Example of outfile.quartiles.exons-introns.definition.sizes.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.quartiles.exons-introns.definition.sizes.png">
</p>
[Figure quartiles exons-introns definition sizes]


- **outfile.density.distribution.[exons|introns].genome.ps** \
These images show some measurements regarding the exon (up image) and intron (down image) density.

    On the left axis, it is shown the distribution of the total or filtered number (between the lower fence and the upper fence) of exons (up image) or introns (down image) according to their position within the protein-coding gene. On the rigth axis, only for the intron's image, the ratio of introns per exons is plotted by taking into account either all exons and introns (standard ratio) or only those exons and introns between the lower fence and the upper fence (weighted ratio).

    The mean density of all exons per protein-coding gene (the standard average density) as well as of those exons between the lower fence and the upper fence (the filtered average density) are also shown (up image). The mean density of all introns per intron-bearing genes (the standard average density) as well as of those introns between the lower fence and the upper fence (the filtered average density) are also shown (down image).

*Example of outfile.density.distribution.exons.genome.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.density.distribution.exons.genome.png">
</p>
[Figure density distribution exons genome]

*Example of outfile.density.distribution.introns.genome.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.density.distribution.introns.genome.png">
</p>
[Figure]

- **total.outfile.[exons|introns|intergenics]_ranges.table.txt** \
This file provides the proportion (%), from the total population number, of the corresponding feature sizes (exons, introns or intergenics) by nucleotide ranges: 1 to 15, 16 to 50, 51 to 100, 101 to 250 nts, etc.

*Example of the total.outfile.exons_ranges.table.species.txt:*

```terminal
SPECIES   (number)   GENOME_MBS   15   50   100   250   500   1000   5000   10000   50000   100000   	500000   1000000   5000000   10000000   50000000   100000000   500000000	  1000000000
volvox_carteri   (108220)   131.16   0.09   6.27   22.78   48.73   14.69   5.07   2.35   0.03   0   0 0   0   0   0   0   0   0   0
```

*Example of the total.outfile.introns_ranges.table.species.txt:*

```terminal
SPECIES   (number)   GENOME_MBS   15   50   100   250   500   1000   5000   10000   50000   100000	500000   1000000   5000000   10000000   50000000   100000000   500000000   1000000000
volvox_carteri   (92822)   131.16   0.02   1.67   7.80   23.69   40.47   21.44   4.91   0   0   0           	0   0   0   0   0   0   0   0
```

- **outfile.ranges.sizes.distribution.genome.ps** \
The image shows the propotion (%) from the total population number of the corresponding feature sizes (exons, introns or intergenics) by nucleotide ranges: 1 to 15, 16 to 50, 51 to 100, 101 to 250 nts, etc.

*Example of outfile.ranges.sizes.distribution.genome.ps:*

<p align="center">
<img src="/volvox_carteri_outfiles/volvox.carteri.archaeplastida.chlorophyta.phytozome.assembly.ranges.sizes.distribution.genome.png">
</p>
[Figure density distribution exons genome]

- **total.outfile.intron_types.table.txt | total.outfile.intron_types.strands.table.txt** \
	Each file provides the estimation about the excess/deficit of the intron length distributions modulo 3 (as suggested by Roy and Penny, 2007, PMID:17617639) for all intron sizes in the genome or by strands. Since introns are not expected to respect the coding frame, intron lengths `3n, 3n+1`, and `3n+2` should appear in similar fractions `p3n ≈ p3n+1 ≈ p3n+2`.  As stated by Roy and Penny (2007), large values of the `''**3n excess**''`, `E3 = p3n − (p3n+1+p3n+2) / 2`,  suggest that a considerable fraction of internal exons could have been incorrectly predicted as introns or that there are several ''intron retention'' events. On the other hand, **`a deficit of 3n introns`**, i.e., `E3 ≪ 0`, may suggest that a considerable fraction of 3n introns –lacking of stop codons– may have been mistaken for exons.

*Example of the total.outfile.intron_types.table.txt:*

```terminal
SPECIES			INTRONS		3n		3n+1	3n+2		EXCESS 3n		EXCESS 3n+1
volvox carteri		92822		0.347	  0.328	0.325		0.020			0.003
```

*Example of the total.outfile.intron_types.strands.table.txt*

```terminal
SPECIES 	STRAND	INTRONS		3n  	3n+1	  3n+2		EXCESS 3n		EXCESS 3n+1
volvox carteri	   +	  45662		0.350	 0.324	 0.327		0.024		-0.003
volvox carteri	   -	  47160		0.344	 0.333	 0.323		0.016		0.010
```

- **total.outfile.[exons|introns|intergenics]_ranges.table.txt** \
These files report the nucleotide composition for the features (exons, introns, intergenics) beloging to particular size ranges.

Example of the total.outfile.exons_ranges_atgc.table.txt

```terminal
SIZE_NTS_RANGE 	      NUMBER	FREQ_GC%	 FREQ_AT%	FREQ_N%
1-15			101		52.67		47.33		0
16-50			6780		57.30		42.70		0
51-100			24648		59.11		40.89		0
101-250			52736		61.30		38.70		0
251-500			15895		62.22		37.78		0
501-1000		5491		64.07		35.93		0
1001-5000		2540		65.77		34.23		0
5001-10000		29		66.76		33.24		0
10001-50000		0		0			0			0
```

- **total.outfile.summary.genome.txt**
- **total.outfile.[exons|introns|intergenics].table.txt**
- **total.outfile.[exons|introns|intergenics].strands.table.txt**
- **total.outfile.[exons|introns|intergenics].chromosomes.table.txt** \
These files report several statistic estimators from the whole population for every genome feature (CDS, exons, introns and intergenic regions) at different levels: the whole genome, by strands and by chromosomes (when available). Meaning of the headers is described in the table of the following section.

__Table for the description of the headers in the output files__


| **VARIABLE**                      | **DESCRIPTION**                                                              |
| --------------------------------- | ---------------------------------------------------------------------------- |
| **MEASUREMENTS FOR GENOME SIZE**  |
| **GENOME\_NTS**                   | Total nucleotide length of the genome (assembled or estimated)               |
| **GENOME\_MBS**                   | Total nulceotide length of the genome (assembled or estimated) as Megabases  |
| **LOG\_GENOME\_SIZE\_NTS**        | Logaritmic value (base 10) of the genome size measured in nucleotides        |
| **LOG\_GENOME\_MBS**              | Logaritmic value (base 10) of the genome size measured in Megabases          |
| **MEASUREMENTS FOR CDS**          |
| **CDS\_NUMBER**                   | Total number of protein-coding genes                                         |
| **CDS\_CONTENT\_MBS**             | Total nucleotide coverage of protein-coding genes in the genome as Megabases |
| **CDS\_AVG\_SIZE**                | Average size of protein-coding genes in the genome                           |
| **CDS\_with\_INTRONS**            | Total number of protein-coding genes harboring introns                       |
| **CDS\_with\_INTRONS\_%**         | Fraction of protein-coding genes harboring introns                           |
| **TOTAL\_NUMBER\_TRANSCRIPTS**    |                                                                              |
| **CDS\_with\_TRANSCRIPTS**        |                                                                              |
| **CDS\_with\_TRANSCRIPTS\_%**     |                                                                              |
| **AVG\_NUMBER\_TRANSCRIPT**       |                                                                              |
| **WEIGHTAVG\_NUMBER\_TRANSCRIPT** |                                                                              |
| **MEASUREMENTS for INTRONS, EXONS and INTERGENICS** |
| **TOTAL\_NUMBER**                                   | Total number of introns/exons in the genome -TRANSTRIPT OR CONTENT?-|
| **AVERAGE\_SIZE**                                   | **Straight average:** average of all intron/exon sizes found in the protein-coding genes of a genome. It depends on the number of data points from the whole sample (_i.e._, gene models with exons/introns), which equally contribute to the final average regardless of which gene they belong to. In the case of intergenic regions, the average is calculated by taking into account all intergenic regions.|
| **LOG\_AVG\_SIZE**                                  | Logaritmic value (base 10) of the intron/exon/intergenics average size measured in nucleotides|
| **FILTERAVG\_SIZE**                                 | Straight average size of the introns/exons/intergenics located between the upper and lower fences of the size distribution|
| **WEIGHTAVG\_SIZE**                                 | **Weighted average size for introns:** average of the mean intron sizes for all protein-coding genes bearing introns. **Weighted average size for exons:** average of the mean exon sizes for all protein-coding genes. It depends on the gene-structure of the genome, and thus, it samples more broadly the data points that contribute –in the case of introns– to the very well known high-skewed length distribution. **Not available for intergenic regions.** |
| **LOG\_WAVG\_SIZE**                                 | Logaritmic value (base 10) of the intron/exon weighted average size measured in nucleotides|
| **WEIGHTAVG\_SIZE\_FILTER**                         | Weighted average size of the introns/exons located between the upper and lower fences of the size distribution |
| **NORMALISEDAVG\_SIZE**                             | Mean size of the introns/exons based on the AVERAGE\_SIZE,  WEIGHTAVG\_SIZE and FILTERAVG\_SIZE|
| **AVG\_VAR**                                        | Variance of the straight AVERAGE\_SIZE|
| **AVG\_SD**                                         | Standard deviation of the AVERAGE\_SIZE|
| **AVG\_EBARBEGIN**                                  | First position of the error bar estimated for the AVERAGE\_SIZE|
| **AVG\_EBAREND**                  | AVG_EBAREND	Final position of the error bar estimated for the AVERAGE_SIZE |
| **AVG_VARFILTER**                 | Variance of the straight FILTERAVG_SIZE |
| **AVG_SDFILTER**                  | Standard deviation of the FILTERAVG_SIZE |
| **AVG_EBARBEGIN_FILTER**          | First position of the error bar estimated for the FILTERAVG_SIZE |
| **AVG_EBAREND_FILTER**            | Final position of the error bar estimated for the FILTERAVG_SIZE |
| **WEIGHTAVG_VAR**                 |	Variance of the WEIGHTAVG_SIZE |
| **WEIGHTAVG_SD**                  |	Standard deviation of the WEIGHTAVG_SIZE |
| **WEIGHTAVG_EBARBEGIN**           |	First position of the error bar estimated for the WEIGHTAVG_SIZE |
| **WEIGHTAVG_EBAREND**             |	Final position of the error bar estimated for the WEIGHTAVG_SIZE |
| **WEIGHTAVG_VARFILTER**           |	Variance of the straight WEIGHTAVG_SIZE_FILTER |
| **WEIGHTAVG_SDFILTER**            |	Standard deviation of the WEIGHTAVG_SIZE_FILTER |
| **WEIGHTAVG_EBARBEGIN_FILTER**    |	First position of the error bar estimated for the WEIGHTAVG_SIZE_FILTER |
| **WEIGHTAVG_EBAREND_FILTER**      |	Final position of the error bar estimated for the WEIGHTAVG_SIZE_FILTER |
| **MIN_SIZE**                      |	Minimum size of the intron/exon/intergenic sequence |
| **MAX_SIZE**                      |	Maximum size of the intron/exon/intergenic sequence |
| **SIZE_Q_DOWN**                   |	0.05 percentile of the intron/exon/intergenic population size |
| **SIZE_Q1**                       |	First quartile of the intron/exon/intergenic population size |
| **SIZE_Q2-MEDIAN**                |	Second quartile or median of the intron/exon/intergenic population size |
| **SIZE_Q3**                       |	Third quartile of the intron/exon/intergenic population size |
| **SIZE_Q_UP**                     |	0.95 percentile of the intron/exon/intergenic population size |
| **SIZE_LOWER_FENCE**              |	Lower fence (for outliers) of the intron/exon/intergenic population size, LF = Q1 - 1.5 (IQR) |
| **SIZE_UPPER_FENCE**              |	Upper fence (for outliers) of the intron/exon/intergenic population size, UF = Q3 + 1.5 (IQR) |
| **FULL_LENGTH_NTS**               |	Total nucleotide length of the intronic/exonic/intergenic sequences |
| **FULL_LENGTH_MBS**               |	Total nucleotide length of the intronic/exonic/intergenic sequences as Megabases |
| **CONTENT_NTS**                   |	Total nucleotide coverage of intronic/exonic/intergenic sequences in the genome |
| **CONTENT_MBS**                   |	Total nucleotide coverage of intronic/exonic/intergenic sequences in the genome as Megabases |
| **CONTENT_GENOME_%**              |	Fraction of the genome size that correspond to intronic/exonic/intergenic sequences |
| **LOG_GENOME_CONTENT**            |	Logaritmic value (base 10) of the genome content for a particular feature (intron, exon, intergenic) measured in Megabases |
| **MEASUREMENTS for EXONS and INTRONS ONLY** |
| **RATIO_INTRONS:EXONS**           |	The ratio of total introns per total exons in a genome |
| **EXON_DEFINITION_%**             |	Porportion (%) of internal exons flanked by at least one short intron (≤ 250 nts) |
| **EXONFILTER_DEFINITION_%**       |	Porportion (%) of internal exons located between the upper and lower fences of the size distribution and that are flanked by at least one short intron (≤ 250 nts) |
| **WEIGHTED_DENSITY**              |	Average of the mean intron/exons numbers for all protein-coding genes |
| **STANDARD_DENSITY**              |	**For introns:** average number of introns in protein-coding genes bearing introns. **For exons:** average number of exons in protein-coding genes |
| **DENSITY_MIN**                   |	Minimum number of introns in protein-coding genes |
| **DENSITY_MAX**                   |	Maximum number of introns in protein-coding genes |
| **DENSITY_VAR**                   |	Variance of the STANDARD_DENSITY |
| **DENSITY_SD**                    |	Standard deviation of the STANDARD_DENSITY |
| **DENSITY_EBAR_BEGIN**            |	First position of the error bar estimated for the STANDARD_DENSITY |
| **DENSITY_EBAR_END**              |	Final position of the error bar estimated for the STANDARD_DENSITY |
| **DENSITY_FILTER**                |	Average number of introns and exons within protein-coding genes located between the upper and lower fences of their corresponding size distribution |
| **DENSITY_FILTER_VAR**            |	Variance of the DENSITY_FILTER |
| **DENSITY_FILTER_SD**             |	Standard deviation of the DENSITY_FILTER |
| **DENSITY_FILTER_EBARBEGIN**      |	First position of the error bar estimated for the DENSITY_FILTER |
| **DENSITY_FILTER_EBAREND**        |	Final position of the error bar estimated for the DENSITY_FILTER |
| **NUCLEOTIDE COMPOSITION** |
| **FREQ_A_%**                      |	Frequency of adenine in the sequence |
| **FREQ_T_%**                      |	Frequency of thymine in the sequence |
| **FREQ_G_%**                      |	Frequency of cytosine in the sequence |
| **FREQ_C_%**                      |	Frequency of guanine in the sequence |
| **FREQ_AT_%**                     |	Frequency of A-T in the sequence |
| **FREQ_GC_%**                     |	Frequency of G-C in the sequence |
| **FREQ_Ns_%**                     |	Frequency of undetermined nucleotides in the sequence labeled as N |
| **FREQ_Xs_%**                     |	Frequency of undetermined nucleotides in the sequence labeled as X |
| **FREQ_Zs_%**                     |	Frequency of undetermined nucleotides in the sequence labeled as any other way rather than N and X |


## Manual

## Trouble shooting

The user migth get some of the following ERROR messages while running the `GenomeContent` pipeline:

- **The option -single \<yes | no\> and its accompanying parameters are MANDATORY. Please specify the necessary parameters and try again!**

The user has forgotten to provide the mode to run the program, on a single genome (-single yes) or on several genomes (-single no).

- **The following: -a \<annotations\> -g \<genomes\> -o \<outfiles\> -d \<database\> -t \<status\> are MANDATORY arguments to run GenomeContent.pl on a single genome. Please specify the necessary parameters and try again!**

The user has failed to provided one or all parameters necessary to run the program in a single mode. A similar ERROR message will occur for the non-single mode.

- **ERROR: The genome sequences are not located in '/directory/genomes'! Please correct that and try again.**

The genome sequence is not located in the name and/or the directory provided for the user. Otherwise, the name of the genome sequence **does not match** with any of the following two options: to the name provided by the user for the **single mode** or to the name of the annotation file in the **non-single mode**.

- **ERROR: File 'example.input.file.genome.fasta' doesn't look like a FASTA format and/or a DNA sequence! Please correct that and try again.**

Is the genome sequence in nucleotides and in FASTA format?  The program stops if the fasta format (>) is missing at the beginning of every sequence and/or when the sequence is not given in nucleotides.

- **ERROR: The name of the sequences in the fasta file DO NOT MATCH to the names of the sequences in the gene annotation file '/directory/annotations'! Please correct that and try again.**

This is an important error. The program relies on matching the sequence names of the genome sequence with the sequence names in the gene annotation file. Usually, the format of both files are congruent from the same database and version of the genome project. Please, take into account that files from different databases (e.g., a gene annotation file from ensembl and a genome sequence from phytozome) or even from different versions of  the same genome projects (e.g., human genome GRCh35 versus GRCh38 in ensembl) migth have either incongruencies in the format or in the mapping of the coordinates into the genome sequence when changes in the contigs, scaffolds or chromosomes have occured along the versions.


## Description

`GenomeContent` was written in Perl to calculate global statistics and sequence-based estimators of genome features in six major steps, as shown in Figure 1. First, the processing of gene annotations focuses on identifying coordinates from protein-coding gene (CDS), while the filtering process focuses on checking the quality of intron annotations. As described in next sections, the coordinates derived from both proceses are taken as the “reference gene sets” for introns, exons and intergenic regions to directly estimate several statistic descriptors, such as size, density and number. Then, the “reference gene sets” are projected onto the genome sequence in both strands, so that the nucleotide contents of each genome-feature are calculated according to the definitions described in a section below. Finally, all statistic descriptors obtained with the program are provided as text files, fasta formats and exploratory figures (as shown in Figure 1). `GenomeContent` runs on an entire genome in few minutes or hours, depending on genome size and the number of annotated genes. `GenomeContent` is available upon request during peer-review, and will be openly available after publication.

![Figure 1]()

__Figure 1.__ A flowchart of the `GenomeContent` pipeline. The genome sequence in fasta format and the gene annotation in standard format are required to run the pipeline. Two major steps take place sequentially to process gene annotation coordinates: identification of protein-coding genes, filtering and checking of protein-coding gene annotations, and calculation of statistic descriptors. The coordinates derived from steps one and two provide “reference sets” for introns, exons and intergenic regions. These sets are then projected onto the genome sequence in both strands and the genome-feature content is calculated according to the equations described below.

#### 1.1. Filtering of gene annotations

The filtering process of `GenomeContent` involves: a) identification of CDS, b) removal of small sizes, c) treatment of isoforms, and d) estimation of systematic errors in CDS. First, only genome coordinates from CDS were extracted, but their corresponding untranslated regions (UTRs) are not included because these are not fully annotated in most genome projects (Davuluri, 2001, PMID:11726928; Stormo, 2000, PMID:10779479). By default, `GenomeContent` exclude introns and exons smaller than 15 nucleotides (nts). Although, small exons and introns (<20 nts) have been reported in several eukaryotes (Gilson, 2006, PMID:16760254; Lim, 2001, PMID: 11572975; Deutsch, 1999, PMID:10454621; Wang, 2006, PMID:16632598; Kupfer, 2004, PMID:15470237; Berget, 1995, PMID:7852296; Volfovsky, 2003, PMID:12799353). Thus, this parameter can be changed by the user. Third, alternative splice variants were kept in the data. To avoid redundant/overestimated data, however, exons with partial or full matching boundaries to exons of other transcripts were overlapped; the same rule was applied to introns. In both cases, their coordinates were joined or replaced accordingly; thus, every exon and intron is only counted once. We call this filtered set of protein-coding gene coordinates for every genome as the “reference gene set”.
Most of the collected gene annotations are based on transcript evidence. Nevertheless, we implemented the approach proposed by Roy and Penny (2007, PMID:17617639) in `GenomeContent` to further estimate systematic errors in CDS annotations by identifying the excess/deficit of the intron-length distributions modulo 3. Since introns are not expected to respect the coding frame, intron lengths 3n, 3n + 1, and 3n + 2 should appear in similar fractions `p3n ≈ p3n + 1 ≈ p3n + 2`. As stated in Roy and Penny (2007), large values of “3n excess”, `E3 = p3n − (p3n + 1 + p3n + 2) /2`, suggest that a considerable fraction of internal exons may have been incorrectly predicted as introns or that there are several “intron retention” events. On the other hand, a deficit of 3 n introns, i.e., `E3 ≪ 0`, suggests that a considerable fraction of 3n introns&#8213;lacking of stop codons&#8213;may have been mistaken for exons.

#### 1.2 Estimation of average sizes and densities of genome features

Several statistical estimators for every genome feature were obtained with GenomeContent. Genome size is defined as the net length of nucleotides and placeholders of all sequences conforming the nuclear genome. The average feature size (e.g., intron size and exon size) of CDSs per genome was calculated in two ways. The straight average size `(Afeature)` is calculated as the total length of all feature sequences (exons or introns) in a genome `(Lfeature)` divided by the total number of all feature sequences (exons or introns, respectively) in a genome `(Nfeature): Afeature = Lfeature /Nfeature`. The straight average depends on the number of data points from the whole sample (i.e., gene models with introns), which equally contribute to the final average regardless of which gene they belong to.

If we now consider a feature to be the average length of the respective feature (exon or intron) within one single gene: a `feature = l feature /n feature`, then the weighted average size, ā f eature , is calculated as the mean of the a f eature values of the respective feature (exons or introns) in a genome, according to:



where n represents the total number of CDSs in a genome when calculating ā exon , or the total number of intron-containing CDSs in a genome when calculating `ā` intron (Gudlaugsdottir, 2007). The weighted average depends on the gene-structure of the genome, and thus it samples more broadly the data points that contribute, in the case of introns, to the well known skewed length distribution.

`GenomeContent` also estimates the abundance of introns within a genome with two different estimates to detect small changes of intron-richness and to buffer dramatic changes among updated genome releases and gene annotations. On the one hand, the abundance of introns across CDSs was estimated as the fraction of intron-containing CDSs from the total number of CDSs (% CDS), which might reflect complete intron loss from CDS structures at the genome level. On the other hand, the abundance of introns and exons within CDSs (“absolute density”) was estimated as the mean number per genome of exons in CDS (exon density), or introns per intron-containing CDSs (intron density), respectively. We employ the “absolute density” because, as observed on Table 1, the average number of introns either per sequence region or from the total number of genes (“normalized density”) depends on both genome size and the number of CDSs, respectively.



> **Table 1. Comparison of intron density estimations: absolute (aID) and normalized (nID).**

For instance, the four species listed in Table 1 exhibit around five introns per CDS when a “normalized” density is estimated from the total number of CDSs. However, the “absolute intron density” clearly shows that, for instance, *X. laevis*, *T. chinensis* and *G. theta* have indeed more introns per CDS on average than the other species, despite of having a smaller fraction of intron-containing CDSs (70.8%, 72.5% and 78.3%, respectively), and lower number of introns and CDS in some cases. Clearly, the bias observed in the “normalized intron density” is produced by larger numbers of total CDS in the genome.

#### 1.3 Estimation of genome contents

GenomeContent also estimates the “feature content” of a given genome, i.e., the proportion of nucleotides of the respective genome features (intron, exon, or intergenic) that contributes to genome size. Since most genome annotations only contain protein-coding regions rather than full transcript models, we count only coding exons and introns delimited by a pair of coding exons. As shown in Figure 1, the program projects the sets of genomic intervals for all exons and introns from coding genes located in the plus strand (set A), the minus strand (set B), and of the isoforms (set C). Since a given nucleotide may be classified differently for different isoforms, we used the following domincance rule in order to obtain a unique classification at the genome level:

> **`Exon > Intron > Intergenic_region`**

It reflects the idea that a genomic position is exonic whenever it appears in a coding exon of at least one transcript. Thus, the exon content of a genome is calculated as the total number of nucleotides in the genome sequence that are classified as coding exon with respect to at least one isoform. Analogously, a position is classified as ’intronic’ if it appears inside the boundaries of annotated coding exons, but it does not overlap with any coding sequence. Thus, intron content was determined as the total number of nucleotides of a genome that were classified as intronic. The CDS content of a genome is calculated as the total number of nucleotides covered by the intronic and exonic positions within coding genes. Finally, the non-coding DNA content was computed analogously as the total number of nucleotides in a genome that are not covered by exonic and intronic positions from coding genes. Genome-feature contents are reported as: total size in Megabases (Mb), fraction (%) from the total genome size, and as genomic coordinates.


---
## Distribution
This information is free. You can redistribute it and/or use it for educational or recreational purpose.

`GenomeContent` is a resource that can be distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. These scripts can be redistributed and/or used and/or modified for your research, educational or recreational purpose.

---
## Credits

## Acknowledgement
