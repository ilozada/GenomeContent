#! /usr/bin/perl -w

##---------------------------------------------------------------------------------------##
##                                                                                       ##
##                                                                                       ##
##                         The GENOME CONTENT program                                    ##
##                                                                                       ##
##  Author	            :  Irma Lozada-Chavez                                        ##
##  Version	            :  0.80                                                      ##
##  Last release date       :  June, 2022                                                ##
##  Please send questions, comments and bug reports to: ilozada@bioinf.uni-leipzig.de    ##
##  Documentation & updates :  https://github.com/ilozada/GenomeContent                  ##
##                                                                                       ##
##                                                                                       ##
##  Purpose =>             A program that estimates global statistics and                ##
##                         sequence-based estimators of genome features                  ##
##                         from gene annotation files in GTF* and GFF* formats           ##
##                                                                                       ##
##                                                                                       ##
##   "GenomeContent" program is written in Perl language to achieve three main goals:    ##
##  (1) To evaluate and correctly parse gene annotations in a systematic                 ##
##      and automatic manner, based on the vocabulary and attribute relationships        ##
##      defined by the Sequence Ontology project.                                        ##
##  (2) To retrieve massive multidimensional datasets from gene annotations,             ##
##      with a primary focus on exons and introns of protein- coding genes,              ##
##      as well as intergenic regions and repeats.                                       ##
##  (3) To calculate global statistics of genome-wide features                           ##
##      with two different approaches: gene-based and genome-based estimators.           ##
##                                                                                       ##
##                                                                                       ##
##   This version of the program is still a draft and being subject to tests.            ##
##   But if you use `GenomeContent` for your research, please cite the following paper:  ##
##                                                                                       ##
##   Lozada-Chávez I, Stadler,P.F., and Prohaska, S.J. (2018)                            ##
##   "Genome-wide features of introns are evolutionary decoupled among themselves        ##
##    and from genome size throughout Eukarya"                                           ##
##   https://www.biorxiv.org/content/early/2018/03/18/283549 (re-submitted)              ##
##   doi: https://doi.org/10.1101/283549                                                 ##
##                                                                                       ##
##                                                                                       ##
##---------------------------------------------------------------------------------------##


use strict;
use Getopt::Long;
use POSIX qw(log10);
use Math::Complex;
use Chart::Gnuplot;
use GD::Graph;


#----- global varibles & usage -----#

my ($single, $annotations, $genomes, $outfiles, $species, $database, $status, $minexonsize, $minintronsize, $help) = ("", "", "", "", "", "", "", "", "", 0, 0, "");
my @files = ();                   my @fasta = ();                    my @sequences = ();
my @y_exons = ();                 my @y_introns = ();                my @y_intergenics = ();
my $correct = "";                 my $error = "";                    my $infile = "";
my $content_introns = "";         my $content_exons = "";            my $content_intergenics = "";
my $sum_exons = 0;                my $sum_introns = 0;               my $sum_intergenics = 0;
my $IntronContentGenome = 0;      my $ExonContentGenome = 0;         my $IntergenicContentGenome = 0; 
my $NfreqGenome_introns = 0;      my $ATfreqGenome_introns = 0;      my $GCfreqGenome_introns = 0;      
my $NfreqGenome_exons = 0;        my $ATfreqGenome_exons = 0;        my $GCfreqGenome_exons = 0;        
my $NfreqGenome_intergenics = 0;  my $ATfreqGenome_intergenics = 0;  my $GCfreqGenome_intergenics = 0;  


#---------------------------------------------------------------------------------------------

my $usage= << "END";
$0: is a program that estimates the genomic content and several statistic descriptors 
    for introns, exons and intergenic regions according to the protein-coding gene annotation 
    and genome sequence provided in standard formats 


  USAGE:   
        one genome:      perl $0 -s <single y> -a <annotations> -g <genomes> -o <outfiles> -d <database> -t <status> -e <species>
	several genomes: perl $0 -s <single n> -a <annotations> -g <genomes> -o <outfiles>


  Example for one genome:       perl $0  -s y  -a /home/user/volvox_carteri.genes.gff  -g /home/user/volvox_carteri.assembly.fasta 
                                         -o /home/user/outputs  -d ensembl  -t assembly  -e <species>

  Example for several genomes:  perl $0  -a /home/user/annotations -g /home/user/genomes -o /home/user/outputs
  

  Options:

  -s <single>	   Run $0 for one single genome project: yes [YyTt]
  		   - Requires the path and name of the annotation file, e.g.: -a /home/user/volvox_carteri.genes.[gff/gtf]
		   - Requires the path and name of the genome file, e.g.: -g /home/user/volvox_carteri.assembly.[fa/faa/fasta]

  		   Run $0 for several genome projects, i.e., not in a single mode: no [NnFf]
  		   - Requires the directory where all annotations are, e.g.: -a /home/user/annotations
                   - Requires the directory where all genomes are, e.g.: -g /home/user/genomes
		   - Requires the directory to write the output files, e.g.: -o /home/user/outputs
		   - The names of the annotation and genome files should have the following format:
		     species_name.database.status[genome/assembly].remaining.extension
		        Example for annotations: volvox_carteri.ensembl.genome.whatever_else.[gff/gtf]
		        Example for genomes:     volvox_carteri.ensembl.genome.whatever_else.[fa/faa/fasta]


  -a <annotations>  SINGLE mode:  path and name of the annotation file in GFF or GTF format, e.g., 
                                  -a /home/user/volvox_carteri.genes.gff
		    SEVERAL mode: directory where all annotations are, e.g.: -a /home/user/annotations

  -g <genome>	    SINGLE mode:  path and name of the genome file in FASTA format, e.g.,
                                  -g /home/user/volvox_carteri.assembly.fasta
                    SEVERAL mode: directory where all genomes are, e.g.: -g /home/user/genomes

  -o <outfiles>	    Directory to write the output files, e.g.: -o /home/user/outputs

  -d <database>	    [REQUIRED] for SINGLE mode: Name of the database where the annotation file was downloaded, 
                    e.g.: ensembl, phytozome, jgi, personal, etc

  -t <status>	    [REQUIRED] for SINGLE mode: Status of the sequenced genome, e.g.:
                    genome  = if complete until chromosome assignments
		    assmbly = if assambled at the conting and/or scaffold level 

  -e <species>      [REQUIRED] for SINGLE mode: species name of the genome project, e.g.:
                    volvox, volvox_carteri, vcn

  -h <help>	    Printing this help message!

END
#------------------------------------------------------------



GetOptions(
	"s|single=s"          => \$single,
	"a|annotations=s"     => \$annotations,
        "g|genomes=s"	      => \$genomes,
	"o|outfiles=s"	      => \$outfiles,
	"e|species=s"	      => \$species,
	"d|database=s"        => \$database,
	"t|status=s"          => \$status,
	"me|minexonsize=i"    => \$minexonsize,
	"mi|minintronsize=i"  => \$minintronsize,
	"h|help=s"	      => \$help,
);

##------------------------------------------------------------



#######################################
#                                     #
#   VALIDATING THE INPUT PARAMETERS   #
#                                     # 
#######################################

if ( $help || !$single) {
    die "\n$usage\nThe option -s <single> and its accompanying parameters are mandatory.\n== Please specify the necessary parameters and try again! ==\n";
}

if ( $single =~ /[YyTt]/ ) {
    ($correct, $error)= &CheckYesSingleOptions();
    if ( !$correct ) { die "\nERROR: $error /"; }

    my @split = (); my $split = ""; my $annotationFile = "";
    $split = $annotations;
    $annotations = "";
    $split =~ s/\/\s+$//gis;
    @split = split "/", $split;
    $annotationFile = pop @split;
    $annotations = join ("/", @split);

    @files = ();
    chdir $annotations;
    $files[0] = $annotationFile;


    @split = (); $split = ""; my $sequenceFile = "";
    $split = $genomes;
    $genomes = "";
    $split =~ s/\/\s+$//gis;
    @split = split "/", $split;
    $sequenceFile = pop @split;
    $genomes = join ("/", @split);

    @fasta = ();
    opendir DIRTWO, "$genomes";
    @fasta = readdir DIRTWO;
    closedir DIRTWO;
}


if ( $single =~ /[NnFf]/ ) {
    ($correct, $error)= &CheckNoSingleOptions();
    if ( !$correct ) { die "\nERROR: $error /"; }

    print "this is NO single option\n";

    @files = ();
    chdir $annotations;
    opendir DIR, "$annotations";
    @files = readdir DIR;
    @files = sort {lc($a) cmp lc($b)} @files;
    closedir DIR;

    @fasta = ();
    opendir DIRTWO, "$genomes";
    @fasta = readdir DIRTWO;
    closedir DIRTWO;
}

sub CheckYesSingleOptions {
    if ( !$annotations || !$genomes || !$database || !$status || !$species ) {
        print STDERR $usage;
	return (0, "-a <annotations> -g <genomes> -o <outfiles> -d <database> -t <status> are MANDATORY arguments to run $0 on a single genome\n== Please specify the necessary parameters and try again! ==\n");
    }
    else { 
        print STDERR "- This is the file for the gene annotations: $annotations\n";
        print STDERR "- This is the file for the genome sequences: $genomes\n";

        if ( !$outfiles ) {
            $outfiles = getcwd();
            print STDERR "- The outfiles will be saved on: $outfiles \n";
        }
	else { $outfiles =~ s/\/$//; print STDERR "- This is the directory where the outfile will be saved: $outfiles\n"; }
        return (1, "");
    }
}

sub CheckNoSingleOptions {
    if ( !$annotations || !$genomes || !$outfiles ) {
        print STDERR $usage;
	return (0, "-a <annotations> -g <genomes> -o <outfiles> are MANDATORY arguments to run $0 on several genomes\n== Please specify the necessary parameters and try again! ==\n");
    }
    else {
        print STDERR "- This is the directory for the gene annotations: $annotations\n";
        print STDERR "- This is the directory for the genome sequences: $genomes\n";

	if ( !$outfiles ) {
            $outfiles = getcwd();
            print STDERR "- The outfiles will be saved on: $outfiles \n";
        }
	else { $outfiles =~ s/\/$//; print STDERR "- This is the directory where the outfile will be saved: $outfiles\n"; }
        return (1, "");
    }
}

if ( !$minexonsize ) { $minexonsize = 15; }
if ( !$minintronsize ) { $minintronsize = 15; }
if ( ! -e $outfiles ) { mkdir $outfiles; }

#------------------------------------------------------------


#system "rm $outfiles/*.fasta $outfiles/*.txt $outfiles/*.ps";

if ( $single=~ /[NnFf]/ ) { &openOutfilesSeveralOption($outfiles); }


#----- main cycle for each genome project -----#

foreach $infile ( @files ) {
    chomp $infile;

    if ( $infile =~ /^.+\.(gtf)$/ || $infile =~ /^.+\.(gff)$/ ) {
        #print "*$infile*\n";


        #----- starting to declare global genome variables -----#

        my %ChromosomeSize           = ();   my %CoordinatesOutfile    = ();   my %CoordinatesForOverlaps = ();   my %CoordinatesByAnnotations = ();
        my %CDSbyChromosomes         = ();   my %CDSbyStrands          = ();   my %ATCGbySpecies          = ();	  my %NoncodingChromosomes = ();
        my %NumberByChromosomes      = ();   my %NumberByStrands       = ();   my %NumberBySpecies        = ();   my %SizesNoncodingChromosomes = ();
        my %LengthByChromosomes      = ();   my %LengthByStrands       = ();   my %LengthBySpecies        = ();   my %GenomeNTS = ();
        my %NTSlengthByChromosomes   = ();   my %NTSlengthByStrands    = ();   my %NTSlengthBySpecies     = ();   my %ExonPositionSizes = ();
        my %StatisticsByChromosomes  = ();   my %StatisticsByStrands   = ();   my %StatisticsBySpecies    = ();   my %ExonPositionStatistics = ();
        my %IntronPositionBySpecies  = ();   my %ExonPositionBySpecies = ();   my %IntronPositionSizes    = ();   my %IntronPositionStatistics = ();
        my %WeightAvgByChromosomes   = ();   my %WeightAvgByStrands    = ();   my %WeightAvgBySpecies     = ();   my %Introns = ();
        my %CDSintronsByChromosomes  = ();   my %CDSintronsByStrands   = ();   my %CDSintronsBySpecies    = ();   my %Exons = ();
        my %IntronClassByChromosomes = ();   my %IntronClassByStrands  = ();   my %IntronClassBySpecies   = ();   my %ATGCexons = ();
        my %ATGCintergenics          = ();   my %Intergenics           = ();   my %ATGCintrons            = ();   my %numberExonsIntronsByCDSgenome = ();
        my %meanByCDSbyChromosomes   = ();   my %meanByCDSbyStrands    = ();   my %meanByCDSbySpecies     = ();

        my $ref_hashByCoordinates    = {};   my $ref_hashByAnnotations = {};   my $ref_hashByATGCgenome   = {};   my $ref_hashByATGCintergenics = {};
        my $ref_hashByFasta          = {};   my $ref_hashByExons       = {};   my $ref_hashByATGCexons    = {};
        my $ref_hashByIntergenics    = {};   my $ref_hashByIntrons     = {};   my $ref_hashByATGCintrons  = {};

        my @coordinatesCDS = ();   my @x_intronup   = ();   my @Exons       = ();   my @cds         = ();
        my @cdsSizeGenome  = ();   my @y_introndown = ();   my @Introns     = ();   my @intergenics = ();
        my @definition_set = ();   my @z_exon       = ();   my @Intergenics = (); 

        my $GenomeFrequencyMbs = 0;       my $lengthSequence = 0;            my $geneSize = 0;             my $Agenome = 0;
        my $TotalCDsSizeGenome = 0;       my $lengthExonsByCDS = 0;          my $exon_intron = 0;          my $Tgenome = 0;
        my $IntronsExonsGenome = 0;       my $numberIntronsByCDS = 0;        my $numberExonsByCDS = 0;     my $Ggenome = 0;
        my $lengthIntronsByCDS = 0;       my $totalIntronSize = 0;           my $meanIntronsByCDS = 0;     my $Cgenome = 0;
        my $CDsSizeGenomeMbs = 0;         my $totalExonSize = 0;             my $meanExonsByCDS = 0;       my $Ngenome = 0;
        my $TotalCDsGenome = 0;           my $GenomeSize = 0;                my $cdsSize = 0;              my $ATgenome = 0;
        my $AvgCDsGenome = 0;             my $intron_up = 0;                 my $intron_down = 0;          my $GCgenome = 0;
        my $errorbarBegin = 0;            my $errorbarEnd = 0;               my $errorbarBeginFilter = 0;  my $errorbarEndFilter = 0;
        my $densityIntronsByExons = 0;    my $intergenicNumber = 0;          my $intronNumber = 0;         my $exonNumber = 0;

        my $firstPositionExon = 0;        my $secondPositionExon = 0;        my $firstPositionIntron = 0;  my $secondPositionIntron = 0;
        my $firstPositionIntergenic = 0;  my $secondPositionIntergenic = 0;  my $intergenicSize = 0;       my $lower_fence = 0;
        my $lengthIntergenic = 0;         my $lengthIntron = 0;              my $lengthExon = 0;           my $upper_fence = 0;
        my $exon_range = 0;               my $avg_filter = 0;                my $median = 0;               my $variance = 0;
        my $stddev = 0;                   my $var_filter = 0;                my $sd_filter = 0;            my $null = 0;
        my $Q1 = 0;                       my $Q3 = 0;                        my $Q_top = 0;                my $Q_down = 0;
        my $IQR = 0;                      my $total_nts = 0;		     my $empty_cds = 0;		   my $seqs = 0;
        my $CDSbySpecies = 0;             my $existing = 0;

        my $atExon = "";   my $atIntron = "";   my $atIntergenic = "";   my $seqann          = "";   my $sequenceFile   = "";
        my $gcExon = "";   my $gcIntron = "";   my $gcIntergenic = "";   my $seqname         = "";   my $annotationFile = "";
        my $aExon  = "";   my $aIntron  = "";   my $aIntergenic  = "";   my $coordinatesCDS  = "";   my $check          = "";
        my $tExon  = "";   my $tIntron  = "";   my $tIntergenic  = "";   my $exon_definition = "";   my $nt             = "";
        my $gExon  = "";   my $gIntron  = "";   my $gIntergenic  = "";   my $multipleIntron  = "";   my $gffout         = "";
        my $cExon  = "";   my $cIntron  = "";   my $cIntergenic  = "";   my $feature         = "";   my $rest           = "";
        my $nExon  = "";   my $nIntron  = "";   my $nIntergenic  = "";   my $strand          = "";   my $id             = "";

        my ($g, $h, $i, $j, $k, $l, $m, $n, $o, $p, $q, $r, $s, $t, $v, $w, $x, $y, $z) = 
           (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        my ($aa, $bb, $cc, $dd, $ee, $ff, $gg, $ii, $jj, $kk, $ll, $mm, $nn, $oo, $pp, $qq, $rr, $ss, $xx, $yy, $zz, $ww) = 
           (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        #----- ending of declare global genome variables -----#


        #----- starting to open individual output files -----#

	if ( $single=~ /[YyTt]/ ) { $annotationFile = $infile; }

	if ( $single=~ /[NnFf]/ ) { 
	    $annotationFile = $infile;
	    ($species, $database, $status, $rest) = ("", "", "", "");
	    $infile =~ s/[\.]/=/gis;
	    ($species, $database, $status, $rest) = split "=", $infile;
        }
	$gffout = "$outfiles/$species.$status";

	open SIZE, ">$gffout.size.txt" || die "Cannot open size.out file\n";
	open ALL, ">$gffout.coordinates.txt" || die "Cannot open coordinates.out file\n";
	open FILTER_EXONS, ">$gffout.filtered.exons.fasta" || die "Cannot open exons.out file\n";
        open FILTER_INTRONS, ">$gffout.filtered.introns.fasta" || die "Cannot open introns.out file\n";

	open COORDINATES_EXONS, ">$gffout.exons.content.txt" || die "Cannot open coordinates.exons.content.txt\n";
	open COORDINATES_INTRONS, ">$gffout.introns.content.txt" || die "Cannot open coordinates.introns.content.txt\n";
	open COORDINATES_INTERGENICS, ">$gffout.intergenics.content.txt" || die "Cannot open coordinates.intergenics.content.txt\n";

	open SIZESEXONS_FILTER, ">$gffout.quartiles.sizes.exons.distribution.genome.txt" || die "Cannot open sizes.exons.distribution.txt\n";
	open SIZESINTRONS_FILTER, ">$gffout.quartiles.sizes.introns.distribution.genome.txt" || die "Cannot open sizes.introns.distribution.txt\n";
	open SIZESINTERGENICS_FILTER, ">$gffout.quartiles.sizes.intergenics.distribution.genome.txt" || die "Cannot open sizes.intergenics.distribution.txt\n";

	open SIZESEXONS_ALL, ">$gffout.all.sizes.exons.distribution.genome.txt" || die "Cannot open sizes.exons.distribution.txt\n";
	open SIZESINTRONS_ALL, ">$gffout.all.sizes.introns.distribution.genome.txt" || die "Cannot open sizes.introns.distribution.txt\n";
	open SIZESINTERGENICS_ALL, ">$gffout.all.sizes.intergenics.distribution.genome.txt" || die "Cannot open sizes.intergenics.distribution.txt\n";
	print SIZESEXONS_FILTER "EXON_SIZE\tFREQUENCY\n"; print SIZESINTRONS_FILTER "INTRON_SIZE\tFREQUENCY\n"; print SIZESINTERGENICS_FILTER "INTERGENIC_SIZE\tFREQUENCY\n";
	print SIZESEXONS_ALL "EXON_SIZE\tFREQUENCY\n"; print SIZESINTRONS_ALL "INTRON_SIZE\tFREQUENCY\n"; print SIZESINTERGENICS_ALL "INTERGENIC_SIZE\tFREQUENCY\n";

	open DEFINITION, ">$gffout.order.intron-exon-intron.distribution.txt" || die "Cannot open order.intron-exon-intron.distribution.txt\n";
	open ORDER_EXONS, ">$gffout.order.exons.distribution.txt" || die "Cannot open order.exons.distribution.sizes.txt\n";
	open ORDER_INTRONS, ">$gffout.order.introns.distribution.txt" || die "Cannot open order.introns.distribution.sizes.txt\n";
	print DEFINITION "SPECIE:ID\tEXON_NUMBER\tEXON_SIZE\tINTRON-UP_SIZE\tINTRON-DOWN_SIZE\n";
        print ORDER_EXONS "EXON_POSITION\tTOTAL_EXONS\tSIZE_AVE\tSIZE_FILTER\tMIN\tMAX\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\n";
        print ORDER_INTRONS "INTRON_POSITION\tTOTAL_INTRONS\tSIZE_AVE\tSIZE_FILTER\tMIN\tMAX\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\n";

	open DEFINITION_FILTER, ">$gffout.quartiles.order.intron-exon-intron.distribution.txt" || die "Cannot open intron-exon-intron.distribution.txt\n";
	open ORDEREXONS_FILTER, ">$gffout.quartiles.order.exons.distribution.txt" || die "Cannot open exons.distribution.sizes.txt\n";
	open ORDERINTRONS_FILTER, ">$gffout.quartiles.order.introns.distribution.txt" || die "Cannot open introns.distribution.sizes.txt\n";
	print DEFINITION_FILTER "SPECIE:ID\tEXON_NUMBER\tEXON_SIZE\tINTRON-UP_SIZE\tINTRON-DOWN_SIZE\n";
        print ORDEREXONS_FILTER "EXON_POSITION\tTOTAL_EXONS\tSIZE_AVE\tMIN\tMAX\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\n";
        print ORDERINTRONS_FILTER "INTRON_POSITION\tTOTAL_INTRONS\tSIZE_AVE\tMIN\tMAX\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\n";

        #----- ending of opening individual output files -----#


        #----- getting annotation coordinates by sequences -----#

	my %CoordinatesByChromosomes = ();
	($ref_hashByCoordinates, $ref_hashByAnnotations) = &getCoordinates("$annotations/$annotationFile", $database);
	%CoordinatesByChromosomes = %{$ref_hashByCoordinates};
	%CoordinatesByAnnotations = %{$ref_hashByAnnotations};


        #----- getting genome sequences -----#

 	if ( $single=~ /[YyTt]/ ) { &openOutfilesSingleOption($outfiles, $species); }

        @sequences = ();
	@sequences = grep { /^$species\..+\.fa[a|sta]*$/ && -f "$genomes/$_" } @fasta;

        if ( scalar @sequences == 0 ) {
	    die " ERROR: The genome sequences are not located in '$genomes'!\n Please correct that and try again.";
	    exit;
        }


	my $FastaByChromosomes = {};
	foreach $sequenceFile ( @sequences ) {
	    #print "$sequenceFile\n";
	    my $ref_hashByFasta = {};
	    $ref_hashByFasta = &cutFastaSequences("$genomes/$sequenceFile");
	    @$FastaByChromosomes{keys %$ref_hashByFasta} = values %$ref_hashByFasta;
	}
        @sequences = (); $ref_hashByFasta = {};
	#for $seqname ( sort keys %$FastaByChromosomes ) { print "$seqname\t$$FastaByChromosomes{$seqname}\n"; }
	$seqs = scalar keys %$FastaByChromosomes;
	print "Total number of sequences in the genome: $seqs\n";


	#----- cheking point for the sequence names in both files: fasta and annotation -----#

	foreach $seqann ( sort keys %CoordinatesByAnnotations ) {
            if ( $$FastaByChromosomes{$seqann} ) { $existing++; }
	}
	if ( $existing == 0 ) {
	    die " ERROR: The name of the sequences in the fasta file DO NOT MATCH the names of the sequences in the gene annotation file '$annotations'!\n Please correct that and try again.";
	    exit;
	}


        #----- first main module for each genome -----#
        $IntronContentGenome = 0;      $ExonContentGenome = 0;         $IntergenicContentGenome = 0; 
        $NfreqGenome_introns = 0;      $ATfreqGenome_introns = 0;      $GCfreqGenome_introns = 0;      
        $NfreqGenome_exons = 0;        $ATfreqGenome_exons = 0;        $GCfreqGenome_exons = 0;        
        $NfreqGenome_intergenics = 0;  $ATfreqGenome_intergenics = 0;  $GCfreqGenome_intergenics = 0;  

	$aa = 0; $dd = 0; $gg = 0; $ll = 0; $mm = 0; $nn = 0; $oo = 0; $pp = 0; $qq = 0; $rr = 0; $ss = 0;
	$feature = "CDS";


	NEXT:
	for $seqname ( sort keys %$FastaByChromosomes ) {
	    my $lengthSequence = ""; my $strand = "";
	    $xx = 0; $yy = 0; $zz = 0; $bb = 0; $ee = 0;
	    $o = 0; $p = 0; $q = 0;

	    ($ref_hashByATGCgenome, $lengthSequence) = &freqATGC(\$$FastaByChromosomes{$seqname});
	    $GenomeNTS{$species}{A}+= ${$ref_hashByATGCgenome}{A};
	    $GenomeNTS{$species}{T}+= ${$ref_hashByATGCgenome}{T};
	    $GenomeNTS{$species}{G}+= ${$ref_hashByATGCgenome}{G};
	    $GenomeNTS{$species}{C}+= ${$ref_hashByATGCgenome}{C};
	    $GenomeNTS{$species}{N}+= ${$ref_hashByATGCgenome}{N};

	    $ChromosomeSize{$species}{$seqname} = $lengthSequence;
	    $GenomeSize+= $lengthSequence;
            #print "$seqname ---> $lengthSequence\n$$FastaByChromosomes{$seqname}\n";
	    print SIZE "$seqname\t$lengthSequence\n";


            #----- counting for sequences with no coordinates for protein-coding genes -----#

	    if ( !defined $CoordinatesByChromosomes{$seqname}{"+"}{$feature} && !defined $CoordinatesByChromosomes{$seqname}{"-"}{$feature} ) {
	        #print "$seqname = $lengthSequence\n";
	        $NoncodingChromosomes{$species}+= $lengthSequence;
		$SizesNoncodingChromosomes{$species}{$seqname}{nts} = $lengthSequence;
		$SizesNoncodingChromosomes{$species}{$seqname}{mbs} = sprintf("%.2f",($lengthSequence / 1000000));

	        # To avoid undefined values
	        $CDSbyChromosomes{$species}{$seqname} = 0;
	        $NumberByChromosomes{$species}{$seqname}{exons} = 0;
	        $NumberByChromosomes{$species}{$seqname}{introns} = 0;
                $NumberByChromosomes{$species}{$seqname}{intergenics} = 0;

		next NEXT;
	    }

            #----- additional module for each strand ( + and - ) -----#

            for $strand ( sort keys %{ $CoordinatesByChromosomes{$seqname} } ) {
	        @intergenics = ();
		$ii = 0; $jj = 0; $kk = 0; $ww = 0; $cc = 0; $ff = 0;
		$r = 0; $s = 0; $t = 0;
		#print "*$strand*\n";

                for $id ( sort keys %{ $CoordinatesByChromosomes{$seqname}{$strand}{$feature} } ) {
		    @coordinatesCDS = (); @Exons = (); @Introns = ();
		    %Exons = (); %Introns = (); %ATGCexons = (); %ATGCintrons = (); %IntronPositionSizes = (); %ExonPositionSizes = ();
		    $ref_hashByExons = {}; $ref_hashByIntrons = {}; $coordinatesCDS = ""; $ref_hashByATGCexons = {}; $ref_hashByATGCintrons = {};
		    $meanExonsByCDS = 0; $meanIntronsByCDS = 0; $lengthExonsByCDS = 0; $lengthIntronsByCDS = 0; $totalExonSize = 0; $totalIntronSize = 0; 
		    $cdsSize = 0; $densityIntronsByExons = 0; $numberExonsByCDS = 0; $numberIntronsByCDS = 0;

		    @coordinatesCDS = keys %{ $CoordinatesByChromosomes{$seqname}{$strand}{$feature}{$id} };
		    $CDSbyChromosomes{$species}{$seqname}++;
		    $CDSbyStrands{$species}{$strand}++;
		    $CDSbySpecies++;
		    #print "$seqname : $lengthSequence *$id* -> $strand = @coordinatesCDS\n";

		    ($totalExonSize, $totalIntronSize, $coordinatesCDS, $ref_hashByExons, $ref_hashByIntrons, $ref_hashByATGCexons, $ref_hashByATGCintrons, $cdsSize) = &getExonsIntrons($minexonsize, $minintronsize, $species, "$id|location=$seqname", $gffout, $strand, \$$FastaByChromosomes{$seqname}, \@coordinatesCDS);

		    if ( $cdsSize > 0 ) { $cdsSizeGenome[$gg++] = $cdsSize; $TotalCDsSizeGenome+= $cdsSize; }
		    $intergenics[$ww++] = $coordinatesCDS;
		    $CoordinatesForOverlaps{$species}{$seqname}{intergenics}[$o++] = $coordinatesCDS;
		    #print "$seqname\t$id\t$strand\t$coordinatesCDS --> $check\n";


		    if ( scalar keys %$ref_hashByExons > 0 ) {
		        %Exons = %{$ref_hashByExons};
		        %ATGCexons = %{$ref_hashByATGCexons};
		        $NumberByChromosomes{$species}{$seqname}{exons}+= scalar(keys %Exons);
		        $NumberByStrands{$species}{exons}{$strand}+= scalar(keys %Exons);
		        $NumberBySpecies{$species}{exons}+= scalar(keys %Exons);
		        $check = scalar(keys %Exons);
			$empty_cds ++;

		        if ( $strand eq "+" ) { @Exons = sort {$a <=> $b} keys %Exons; }
		        if ( $strand eq "-" ) { @Exons = sort {$b <=> $a} keys %Exons; }

		        foreach $exonNumber ( @Exons ) {
                            ($firstPositionExon, $secondPositionExon, $lengthExon, $atExon, $gcExon, $aExon, $tExon, $gExon, $cExon, $nExon, $geneSize) = (0, 0, 0, "", "", "", "", "", "", "", 0);

		            ($firstPositionExon, $secondPositionExon, $lengthExon, $aExon, $tExon, $gExon, $cExon, $atExon, $gcExon, $nExon, $geneSize) = split ":", $Exons{$exonNumber};
		            #print "EXONS >$sp:$id|CDS($geneSize bps) |location=$seqname:$firstPositionExon-$secondPositionExon:$strand:CDS_exon_$exonNumber |organism=$species|cds=$geneSize(nts):exonic=$totalExonSize(nts)|exon_length=$lengthExon|A=$aExon:T=$tExon:G=$gExon:C=$cExon:AT=$atExon:GC=$gcExon:N=$nExon\n";

		            $lengthExon =~ s/\(.+\)$//; $atExon =~ s/\(.+\)$//; $gcExon =~ s/\(.+\)$//; $nExon =~ s/\(.+\)$//;

			    $ExonPositionBySpecies{$species}{$exonNumber}{summatory} += $lengthExon;
			    $ExonPositionBySpecies{$species}{$exonNumber}{number} ++;
			    $ExonPositionSizes{$species}{$exonNumber} = $lengthExon;
			    push @{ $ExonPositionStatistics{$species}{$exonNumber} }, "$lengthExon";

			    $CoordinatesOutfile{$seqname}{exons}{$strand}[$r++] = "$seqname\t$id\tCDS_exon_$exonNumber\t$firstPositionExon\t$secondPositionExon\t$strand=$lengthExon";
			    $CoordinatesForOverlaps{$species}{$seqname}{exons}[$p++] = "$firstPositionExon:$secondPositionExon";

			    $LengthByChromosomes{$species}{$seqname}{exons}+= $lengthExon;
			    $LengthByStrands{$species}{exons}{$strand}+= $lengthExon;
			    $LengthBySpecies{$species}{exons}+= $lengthExon;

			    foreach $nt ( sort keys %ATGCexons ) {
		                $NTSlengthByChromosomes{$species}{$seqname}{exons}{$nt}+= $ATGCexons{$nt};
			        $NTSlengthByStrands{$species}{exons}{$strand}{$nt}+= $ATGCexons{$nt};
			        $NTSlengthBySpecies{$species}{exons}{$nt}+= $ATGCexons{$nt};
			    }
			    $StatisticsByChromosomes{$species}{$seqname}{exons}[$xx++] = $lengthExon;
			    $StatisticsByStrands{$species}{exons}{$strand}[$ii++] = $lengthExon;
			    $StatisticsBySpecies{$species}{exons}[$ll++] = $lengthExon;
			    $ATCGbySpecies{$species}{exons}[$oo++] = "$lengthExon=$atExon:$gcExon:$nExon";
			    $lengthExonsByCDS+= $lengthExon;						## net length L_ex of all exons in a CDS
		        }
			$numberExonsByCDS = scalar @Exons;
			$numberExonsIntronsByCDSgenome{$species}{exons}[$rr++] = $numberExonsByCDS;
		        $meanExonsByCDS = $lengthExonsByCDS / $numberExonsByCDS;
		        $meanByCDSbyChromosomes{$species}{$seqname}{exons}+= $meanExonsByCDS;
			$meanByCDSbyStrands{$species}{exons}{$strand}+= $meanExonsByCDS;
			$meanByCDSbySpecies{$species}{exons}+= $meanExonsByCDS;
			$WeightAvgBySpecies{$species}{exons}[$dd++] = $meanExonsByCDS;
			$WeightAvgByChromosomes{$species}{$seqname}{exons}[$ee++] = $meanExonsByCDS;
			$WeightAvgByStrands{$species}{exons}{$strand}[$ff++] = $meanExonsByCDS;
		    }

		    if ( scalar keys %$ref_hashByIntrons > 0 ) {
	                %Introns = %{$ref_hashByIntrons};
			%ATGCintrons = %{$ref_hashByATGCintrons};
		        $NumberByChromosomes{$species}{$seqname}{introns}+= scalar(keys %Introns);
			$NumberByStrands{$species}{introns}{$strand}+= scalar(keys %Introns);
			$NumberBySpecies{$species}{introns}+= scalar(keys %Introns);
			$NumberBySpecies{$species}{CDSintrons}+= $numberExonsByCDS;

		        if ( $strand eq "+" ) { @Introns = sort {$a <=> $b} keys %Introns; }
		        if ( $strand eq "-" ) { @Introns = sort {$b <=> $a} keys %Introns; }

		        foreach $intronNumber ( @Introns ) {
			    ($firstPositionIntron, $secondPositionIntron, $lengthIntron, $atIntron, $gcIntron, $aIntron, $tIntron, $gIntron, $cIntron, $nIntron, $geneSize, $multipleIntron) = (0, 0, 0, "", "", "", "", "", "", "", 0, "");

			    ($firstPositionIntron, $secondPositionIntron, $lengthIntron, $aIntron, $tIntron, $gIntron, $cIntron, $atIntron, $gcIntron, $nIntron, $geneSize, $multipleIntron) = split ":", $Introns{$intronNumber};
			    #print "INTRONS >$sp:$id|CDS($geneSize bps) |location=$seqname:$firstPositionIntron-$secondPositionIntron:$strand:CDS_intron_$intronNumber:$multipleIntron |organism=$species|cds=$geneSize(nts):intronic=$totalIntronSize(nts)|intron_length=$lengthIntron|A=$aIntron:T=$tIntron:G=$gIntron:C=$cIntron:AT=$atIntron:GC=$gcIntron:N=$nIntron\n";

			    $lengthIntron =~ s/\(.+\)$//; $atIntron =~ s/\(.+\)$//; $gcIntron =~ s/\(.+\)$//; $nIntron =~ s/\(.+\)$//;

			    if ( $lengthIntron eq "" ) { $lengthIntron = 0; }

			    $IntronPositionBySpecies{$species}{$intronNumber}{summatory} += $lengthIntron;
			    $IntronPositionBySpecies{$species}{$intronNumber}{number} ++;
			    $IntronPositionSizes{$species}{$intronNumber} = $lengthIntron;
			    push @{ $IntronPositionStatistics{$species}{$intronNumber} }, "$lengthIntron";

			    $CoordinatesOutfile{$seqname}{introns}{$strand}[$s++] = "$seqname\t$id\tCDS_intron_$intronNumber\t$firstPositionIntron\t$secondPositionIntron\t$strand=$lengthIntron";
			    $CoordinatesForOverlaps{$species}{$seqname}{introns}[$q++] = "$firstPositionIntron:$secondPositionIntron";

			    $LengthByChromosomes{$species}{$seqname}{introns}+= $lengthIntron;
			    $LengthByStrands{$species}{introns}{$strand}+= $lengthIntron;
			    $LengthBySpecies{$species}{introns}+= $lengthIntron;
			    
			    foreach $nt ( sort keys %ATGCintrons ) {
		                $NTSlengthByChromosomes{$species}{$seqname}{introns}{$nt}+= $ATGCintrons{$nt};
			        $NTSlengthByStrands{$species}{introns}{$strand}{$nt}+= $ATGCintrons{$nt};
			        $NTSlengthBySpecies{$species}{introns}{$nt}+= $ATGCintrons{$nt};
			    }
			    $StatisticsByChromosomes{$species}{$seqname}{introns}[$yy++] = $lengthIntron;
			    $StatisticsByStrands{$species}{introns}{$strand}[$jj++] = $lengthIntron;
			    $StatisticsBySpecies{$species}{introns}[$mm++] = $lengthIntron;
			    $ATCGbySpecies{$species}{introns}[$pp++] = "$lengthIntron=$atIntron:$gcIntron:$nIntron";

			    $lengthIntronsByCDS+= $lengthIntron;
			    $IntronClassByChromosomes{$species}{$seqname}{$multipleIntron}++;
			    $IntronClassByStrands{$species}{$strand}{$multipleIntron}++;
			    $IntronClassBySpecies{$species}{$multipleIntron}++;
		        }
			$CDSintronsByChromosomes{$species}{$seqname}++;
			$CDSintronsByStrands{$species}{$strand}++;      ### number of CDS with introns by strand
			$CDSintronsBySpecies{$species}++;
			$numberIntronsByCDS = scalar @Introns;
			$numberExonsIntronsByCDSgenome{$species}{introns}[$ss++] = $numberIntronsByCDS;
			$meanIntronsByCDS = $lengthIntronsByCDS / $numberIntronsByCDS;			    
		    	$meanByCDSbyChromosomes{$species}{$seqname}{introns}+= $meanIntronsByCDS;
			$meanByCDSbyStrands{$species}{introns}{$strand}+= $meanIntronsByCDS;
			$meanByCDSbySpecies{$species}{introns}+= $meanIntronsByCDS;
			$WeightAvgBySpecies{$species}{introns}[$aa++] = $meanIntronsByCDS;
			$WeightAvgByChromosomes{$species}{$seqname}{introns}[$bb++] = $meanIntronsByCDS;
			$WeightAvgByStrands{$species}{introns}{$strand}[$cc++] = $meanIntronsByCDS;
			$densityIntronsByExons = $numberIntronsByCDS / $numberExonsByCDS;
			$meanByCDSbySpecies{$species}{density}+= $densityIntronsByExons;

			foreach $exonNumber ( @Exons ) {
			    $intron_up = 0; $intron_down = 0; $exon_range = 0; my $intron_upsize = 0; my $intron_downsize = 0;
			    $intron_up = $exonNumber - 1;
			    $intron_down = $exonNumber;
			    $IntronsExonsGenome ++;
			    $intron_upsize = $IntronPositionSizes{$species}{$intron_up} || 0;
			    $intron_downsize = $IntronPositionSizes{$species}{$intron_down} || 0; 

			    #print "exon_number: $exonNumber = $ExonPositionBySpecies{$species}{$exonNumber}{size}*\t";

			    if ( $intron_upsize > 0 && $intron_upsize <= 250 || $intron_downsize > 0 && $intron_downsize <= 250 ) { $exon_intron++; }
			    if ( $exonNumber >= 2 && $exonNumber < $Exons[-1] && $intron_upsize > 0 && $intron_downsize > 0 ) {
				    push(@definition_set, ["$id\t$exonNumber\t$ExonPositionSizes{$species}{$exonNumber}", $IntronPositionSizes{$species}{$intron_up}, $IntronPositionSizes{$species}{$intron_down}]);
				    print DEFINITION "$id\t$exonNumber\t$ExonPositionSizes{$species}{$exonNumber}\t$IntronPositionSizes{$species}{$intron_up}\t$IntronPositionSizes{$species}{$intron_down}\n";
			    }
			}
		    }
	        }

		%Intergenics = (); %ATGCintergenics = (); @Intergenics = (); $ref_hashByIntergenics = {}; $ref_hashByATGCintergenics = {};
		if ( scalar @intergenics > 0) {
		    @intergenics = sort {$a cmp $b} @intergenics;
		    #print "Coordenates of CDS for INTERGENICS: @intergenics\n";
	            ($ref_hashByIntergenics, $ref_hashByATGCintergenics) = &getIntergenicRegions($species, "$id|location=$seqname", $gffout, $strand, $lengthSequence, \$$FastaByChromosomes{$seqname}, \@intergenics);

		    if ( scalar keys %$ref_hashByIntergenics > 0 ) {
		        %Intergenics = %{$ref_hashByIntergenics};
		        %ATGCintergenics = %{$ref_hashByATGCintergenics};

		        $NumberByChromosomes{$species}{$seqname}{intergenics}+= scalar(keys %Intergenics);
		        $NumberByStrands{$species}{intergenics}{$strand}+= scalar(keys %Intergenics);
		        $NumberBySpecies{$species}{intergenics}+= scalar(keys %Intergenics);

		        if ( $strand eq "+" ) { @Intergenics = sort {$a <=> $b} keys %Intergenics; }
		        if ( $strand eq "-" ) { @Intergenics = sort {$b <=> $a} keys %Intergenics; }

		        foreach $intergenicNumber ( @Intergenics ) {
		            ($firstPositionIntergenic, $secondPositionIntergenic, $lengthIntergenic, $atIntergenic, $gcIntergenic, $aIntergenic, $tIntergenic, $gIntergenic, $cIntergenic, $nIntergenic, $intergenicSize) = (0, 0, 0, "", "", "", "", "", "", "", 0);

		            ($firstPositionIntergenic, $secondPositionIntergenic, $lengthIntergenic, $aIntergenic, $tIntergenic, $gIntergenic, $cIntergenic, $atIntergenic, $gcIntergenic, $nIntergenic, $intergenicSize) = split ":", $Intergenics{$intergenicNumber};
		            #print "INTERGENICS >$sp:$seqname:$strand:$firstPositionIntergenic-$secondPositionIntergenic:intergenic_$intergenicNumber($intergenicSize nts) |organism=$species|length=$lengthIntergenic|A=$aIntergenic:T=$tIntergenic:G=$gIntergenic:C=$cIntergenic:AT=$atIntergenic:GC=$gcIntergenic:N=$nIntergenic\n";

		            $lengthIntergenic =~ s/\(.+\)$//; $atIntergenic =~ s/\(.+\)$//; $gcIntergenic =~ s/\(.+\)$//; $nIntergenic =~ s/\(.+\)$//;
		            $CoordinatesOutfile{$seqname}{intergenics}{$strand}[$t++] = "$seqname\tintergenic_$intergenicNumber\t$firstPositionIntergenic\t$secondPositionIntergenic\t$strand=$lengthIntergenic";     

		            $LengthByChromosomes{$species}{$seqname}{intergenics}+= $lengthIntergenic;
		            $LengthByStrands{$species}{intergenics}{$strand}+= $lengthIntergenic;
		            $LengthBySpecies{$species}{intergenics}+= $lengthIntergenic;
		    
		            foreach $nt ( sort keys %ATGCintergenics ) {
		                $NTSlengthByChromosomes{$species}{$seqname}{intergenics}{$nt}+= $ATGCintergenics{$nt};
			        $NTSlengthByStrands{$species}{intergenics}{$strand}{$nt}+= $ATGCintergenics{$nt};
			        $NTSlengthBySpecies{$species}{intergenics}{$nt}+= $ATGCintergenics{$nt};
		            }
		            $StatisticsByChromosomes{$species}{$seqname}{intergenics}[$zz++] = $lengthIntergenic;
		            $StatisticsByStrands{$species}{intergenics}{$strand}[$kk++] = $lengthIntergenic;
		            $StatisticsBySpecies{$species}{intergenics}[$nn++] = $lengthIntergenic;
		            $ATCGbySpecies{$species}{intergenics}[$qq++] = "$lengthIntergenic=$atIntergenic:$gcIntergenic:$nIntergenic";
		        }
		    }
		}
	    }
            unless ( defined $NumberByChromosomes{$species}{$seqname}{introns} ) { $NumberByChromosomes{$species}{$seqname}{introns} = 0; }
        }
	$FastaByChromosomes = {}; @cds = ();
	$GenomeFrequencyMbs = 0; $TotalCDsGenome = 0; $AvgCDsGenome = 0; $CDsSizeGenomeMbs = 0;
	$errorbarBegin = 0; $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;
	$Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0;
        $median = 0; $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $null = 0;
        $Agenome = 0; $Tgenome = 0; $Ggenome = 0; $Cgenome = 0; $Ngenome = 0; $ATgenome = 0; $GCgenome = 0;
	

	$GenomeFrequencyMbs = sprintf("%.2f",($GenomeSize / 1000000));
	$total_nts = $GenomeNTS{$species}{A} + $GenomeNTS{$species}{T} + $GenomeNTS{$species}{G} + $GenomeNTS{$species}{C} + $GenomeNTS{$species}{N};
	#print "AA: $GenomeNTS{$species}{A} -> NTS: $total_nts\n";
	$Agenome = sprintf("%.2f",(($GenomeNTS{$species}{A} * 100) / $total_nts));
	$Tgenome = sprintf("%.2f",(($GenomeNTS{$species}{T} * 100) / $total_nts));
	$Ggenome = sprintf("%.2f",(($GenomeNTS{$species}{G} * 100) / $total_nts));
	$Cgenome = sprintf("%.2f",(($GenomeNTS{$species}{C} * 100) / $total_nts));
	$Ngenome = sprintf("%.2f",(($GenomeNTS{$species}{N} * 100) / $total_nts));
	$ATgenome = sprintf("%.2f",((($GenomeNTS{$species}{A} + $GenomeNTS{$species}{T}) * 100) / $total_nts));
	$GCgenome = sprintf("%.2f",((($GenomeNTS{$species}{G} + $GenomeNTS{$species}{C}) * 100) / $total_nts));
	@cds = sort {$a <=> $b} @cdsSizeGenome;
	$TotalCDsGenome = scalar @cds;
	$AvgCDsGenome = sprintf("%.2f",($TotalCDsSizeGenome / $TotalCDsGenome));
	$CDsSizeGenomeMbs = sprintf("%.2f",($TotalCDsSizeGenome / 1000000));
	#print "CDS SIZES: @cds\n";
	print "Total number of CDS after filtering: $empty_cds\n";

	($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@cds);
	($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null) = &calculateStdDev($AvgCDsGenome, $lower_fence, $upper_fence, \@cds);
	$errorbarBegin = sprintf("%.3f",($AvgCDsGenome - $stddev)); $errorbarEnd = sprintf("%.3f",($AvgCDsGenome + $stddev));
	$errorbarBeginFilter = sprintf("%.3f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.3f",($avg_filter + $sd_filter));
	$exon_definition = "$exon_intron:$IntronsExonsGenome";
	$total_nts = ""; $lengthSequence = "";

	print "**genome -> summary statistics finished**\n";
	print SPECIES "$species\t$GenomeFrequencyMbs\t$CDSbySpecies\t$CDsSizeGenomeMbs\t$AvgCDsGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$Agenome\t$Tgenome\t$Ggenome\t$Cgenome\t$Ngenome\t$ATgenome\t$GCgenome\n";
	#print "SPECIES: $species\t$GenomeFrequencyMbs\t$CDSbySpecies\t$CDsSizeGenomeMbs\t$AvgCDsGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$Agenome\t$Tgenome\t$Ggenome\t$Cgenome\t$Ngenome\t$ATgenome\t$GCgenome\n";

	close EXONS;        close INTRONS;        close INTERGENICS;
	close FILTER_EXONS; close FILTER_INTRONS; close SIZE;


	&getStatisticsBySpecies("exons", $species, $status, $outfiles, $single, $gffout, $GenomeSize, $CDSbySpecies, 
				\%CoordinatesForOverlaps, \%CDSbyChromosomes, \%CDSbyStrands, 
				\%NumberBySpecies, \%NumberByChromosomes, \%NumberByStrands, 
				\%LengthBySpecies, \%LengthByChromosomes, \%LengthByStrands, 
				\%NTSlengthBySpecies, \%NTSlengthByChromosomes, \%NTSlengthByStrands, 
				\%StatisticsBySpecies, \%StatisticsByChromosomes, \%StatisticsByStrands,
				\%ATCGbySpecies, \%ChromosomeSize, \%CoordinatesOutfile, 
				\%meanByCDSbySpecies, \%meanByCDSbyChromosomes, \%meanByCDSbyStrands, 
				\%WeightAvgBySpecies, \%WeightAvgByChromosomes, \%WeightAvgByStrands, \%numberExonsIntronsByCDSgenome);

	&getStatisticsBySpecies("introns", $species, $status, $outfiles, $single, $gffout, $GenomeSize, $CDSbySpecies, 
	                        \%CoordinatesForOverlaps, \%CDSbyChromosomes, \%CDSbyStrands, 
				\%NumberBySpecies, \%NumberByChromosomes, \%NumberByStrands, 
				\%LengthBySpecies, \%LengthByChromosomes, \%LengthByStrands, 
				\%NTSlengthBySpecies, \%NTSlengthByChromosomes, \%NTSlengthByStrands, 
				\%StatisticsBySpecies, \%StatisticsByChromosomes, \%StatisticsByStrands,
				\%ATCGbySpecies, \%ChromosomeSize, \%CoordinatesOutfile, 
				\%meanByCDSbySpecies, \%meanByCDSbyChromosomes, \%meanByCDSbyStrands, 
				\%WeightAvgBySpecies, \%WeightAvgByChromosomes, \%WeightAvgByStrands, \%numberExonsIntronsByCDSgenome, 
				\%ExonPositionBySpecies, \%ExonPositionStatistics, \%IntronPositionBySpecies, \%IntronPositionStatistics,
				\%CDSintronsBySpecies, \%CDSintronsByChromosomes, \%CDSintronsByStrands, \@definition_set, $exon_definition);

	&getStatisticsBySpecies("intergenics", $species, $status, $outfiles, $single, $gffout, $GenomeSize, $CDSbySpecies, 
				\%CoordinatesForOverlaps, \%CDSbyChromosomes, \%CDSbyStrands, 
				\%NumberBySpecies, \%NumberByChromosomes, \%NumberByStrands, 
				\%LengthBySpecies, \%LengthByChromosomes, \%LengthByStrands, 
				\%NTSlengthBySpecies, \%NTSlengthByChromosomes, \%NTSlengthByStrands, 
				\%StatisticsBySpecies, \%StatisticsByChromosomes, \%StatisticsByStrands, 
				\%ATCGbySpecies, \%ChromosomeSize, \%CoordinatesOutfile, 
				\%NoncodingChromosomes, \%SizesNoncodingChromosomes);

	&getIntronModulo3Distributions($species, $outfiles, $single, 
			              \%NumberBySpecies, \%NumberByChromosomes, \%NumberByStrands, 
			              \%IntronClassBySpecies, \%IntronClassByChromosomes, \%IntronClassByStrands);

	close ALL;
	close DEFINITION;        close ORDER_EXONS;         close ORDER_INTRONS;
	close DEFINITION_FILTER; close ORDEREXONS_FILTER;   close ORDERINTRONS_FILTER;
	close SIZESEXONS_ALL;    close SIZESINTRONS_ALL;    close SIZESINTERGENICS_ALL;
	close SIZESEXONS_FILTER; close SIZESINTRONS_FILTER; close SIZESINTERGENICS_FILTER;
	close COORDINATES_EXONS; close COORDINATES_INTRONS; close COORDINATES_INTERGENICS;
    }
}

close SPECIES;			   close STATSTWO_STRAND;	   close STATSTWO_GENOME;
close STATSEXONS_CHROMOSOME;       close STATSEXONS_STRANDS;       close STATSEXONS_GENOME;
close STATSINTRONS_CHROMOSOME;     close STATSINTRONS_STRANDS;     close STATSINTRONS_GENOME;
close STATSINTERGENICS_CHROMOSOME; close STATSINTERGENICS_STRANDS; close STATSINTERGENICS_GENOME;
close EXONSRANGE_GENOME;           close INTRONSRANGE_GENOME;      close INTERGENICSRANGE_GENOME;
close EXONSATGC_GENOME;            close INTRONSATGC_GENOME;       close INTERGENICSATGC_GENOME;





### --------------------------------------- ###
###  *  *  *  S U B R U T I N E S  *  *  *  ###
### --------------------------------------- ###



###########################################################################################################
#                                                                                                         #
#   IDENTIFYING PROTEIN-CODING GENE FEATURES: INTRONS AND EXONS                                           #
#   The goal is to identify full protein-coding genes rather than individual isoforms or transcripts      #
#                                                                                                         #
#   1. Keep only CDS identifiers (usually depicted in the third column of the gtf/gff annotation file)    #
#   2. Only take into account CDS >= 15 nts                                                               #
#   3. Eliminate all unnecessary legends and isoform annotations surrounding the protein-coding gene ID   #
#                                                                                                         # 
###########################################################################################################


sub getCoordinates {
    # INPUT: ($ref_hashByCoordinates, $ref_hashByAnnotations) = &getCoordinates($annotationFile, $database);

    my $db = ""; my $gtffFile = ""; my $line = $_; my %hashByCoordinates = (); my %hashByCoordinatesCheck = (); my %cds = (); my @ids = ();
    my ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others, $id, $cdsLength, $cds) = ("", "", "", 0, 0, 0, "", "", "", "", "", 0, 0);

    $gtffFile = shift; $db = shift;
    open GTFF, "$gtffFile" || die "Cannot open $gtffFile";
    print "$gtffFile\n";


    if ( $db ne "ensembl" ) { 
        while ( <GTFF> ) {
            $line = $_;
	    chomp $line; $line =~ s/\s+$//; $line =~ s/\t+$//;

	    if ( $line !~ /^\#+/ ) {
	        @ids = ();
	        ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others, $id, $cdsLength) = ("", "", "", 0, 0, 0, "", "", "", "", "", 0);
	        ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others) = split "\t", $line;
	        $fastaname =~ s/\s+.+$//;

	        if ( $attributes ne "" && $CDSfeatures eq "CDS" && $source ne "transposable_element" && $source ne "pseudogene" && $source ne "manual_curation_te" && $source ne "tRNA_pseudogene" && $source ne "noncoding_exon" && $source ne "mobile_element" && $source ne "misc_feature" && $source ne "gap" && $source ne "chromosome" && $source !~ /^.+RNA.*/ ) {

                    if ( $firstPosition == 0 ) { $firstPosition = 1;  }
	            $cdsLength = ($secondPosition - $firstPosition) + 1;

		    if ( $cdsLength >= 15 ) {
	                @ids = split ";", $attributes;
	                $id = $ids[0];
	                $id =~ s/^\s+//; $id =~ s/\s+$//;
	                if ( $id =~ /^ID=PAC:/ )						    { $id =~ s/^ID=PAC:/PACid:/; $id =~ s/\.+.+$//; }
	                if ( $id =~ /^ID=CDS:/ && $id =~ /\.\d+$/ )				    { $id =~ s/^ID=CDS://; $id =~ s/\.\d+$//; }
	                if ( $id =~ /^ID=CDS:/ && $id =~ /_\d+$/ )				    { $id =~ s/^ID=CDS://; $id =~ s/_\d+$//; }
	                if ( $id =~ /^ID="CDS:/ && $id =~ /_\d+"$/ )			    	    { $id =~ s/^ID="CDS://; $id =~ s/_\d+"$//; }
		        $id =~ s/^ID=CDS://;
		        $id =~ s/^ID=CDS:PAC://;
	                if ( $id =~ /^ID=/ && $id =~ /-RA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-RA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /-PA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-PA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /-TA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-TA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /(:exon:|:cds:)\d+$/ )			    { $id =~ s/^ID=//; $id =~ s/:exon:\d+$//; $id =~ s/:cds:\d+$//; }
	                if ( $id =~ /^CDS\s+/ )						    	    { $id =~ s/^CDS\s+//; $id =~ s/\s+$//; }
	                if ( $id =~ /^ID=exon_/ )						    { $id =~ s/^ID=exon_//; $id =~ s/-\d+$//; }
	            	$id =~ s/^ID=cds_of_//;
		    	$id =~ s/^ID=gene://; $id =~ s/^ID=mRNA://;
	            	if ( $id =~ /^ID=jCVI\|cds_/ || $id =~ /^ID=gb\|cds_/ || $id =~ /^ID=cds_/ ){ $id =~ s/^ID=jCVI\|cds_//; $id =~ s/^ID=gb\|cds_//; $id =~ s/^ID=cds_//; $id =~ s/-\d+$//; }
	            	if ( $id =~ /^(ID=Trichomonas_vaginalis_Carlton\|cds_)/ )		    { $id =~ s/$1//; $id =~ s/-\d+$//; }
	            	if ( $id =~ /_exon\d+$/ )						    { $id =~ s/_exon\d+$//; $id =~ s/^ID=//; }
	            	if ( $id =~ /_exon_\d+$/ )						    { $id =~ s/_exon_\d+$//; $id =~ s/^ID=//; }
	            	if ( $id =~ /^Parent=/ && $id =~ /\.mRNA\d+$/ )				    { $id =~ s/^Parent=//; $id =~ s/\.mRNA\d+$//; }
	            	if ( $id =~ /^Parent=/ && $id =~ /-RA.*$/ )				    { $id =~ s/^Parent=//; $id =~ s/-RA$//; }
		    	if ( $id =~ /^Parent=/ && $id =~ /-PA.*$/ )			            { $id =~ s/^Parent=//; $id =~ s/-PA$//; }
		    	if ( $id =~ /^Parent=/ && $id =~ /-TA.*$/ )		       	            { $id =~ s/^Parent=//; $id =~ s/-TA$//; }
	            	if ( $id =~ /^GenePrediction\s+/ && $id =~ /-TA.*$/ )		   	    { $id =~ s/^GenePrediction\s+//; $id =~ s/-TA$//; }

	            	$id =~ s/^Name=//; $id =~ s/^name "//; $id =~ s/^gene_id "//; $id =~ s/"+//gis; $id =~ s/^mRNA\s+//gis; $id =~ s/^mRNA://;
	            	$id =~ s/^Name=CDS;Parent=//; $id =~ s/^Parent=//; $id =~ s/^ID=//; $id =~ s/^\|cds_//; $id =~ s/\s+.+$//;
		    	$id =~ s/(-RA|-PA):(exon|cds):\d+$//; $id =~ s/-R[A-Z]$//; $id =~ s/-P[A-Z]$//; $id =~ s/-T[A-Z]$//;

		    	$id =~ s/\.t\d+$//; $id =~ s/\.t\d+\.cds.*$//; $id =~ s/\.t\d+\.CDS.*$//;
		    	$id =~ s/\.mRNA\d+$//; $id =~ s/-mRNA-\d+$//;
		    	$id =~ s/_cds_\d+$//; $id =~ s/\.CDS\d+$//; $id =~ s/:cds$//;
		    	#$id =~ s/\.\d+$//;

	            	#print "IN SUB NON-ENSEMBL = $fastaname = $id -> $CDSfeatures --> $plusORminus:$firstPosition:$secondPosition ---> $cdsLength\n";
	            	$hashByCoordinates{$fastaname}{$plusORminus}{$CDSfeatures}{$id}{"$firstPosition:$secondPosition"} = "$firstPosition:$secondPosition";
		    	$hashByCoordinatesCheck{$fastaname}++;
		    	$cds{$id}++;
		    }
                }
	    }
	}
	close GTFF;
    }

    if ( $db eq "ensembl" ) { 
	while ( <GTFF> ) {
	    $line = $_;
	    chomp $line; $line =~ s/\s+$//; $line =~ s/\t+$//;

	    if ( $line !~ /^\#+/ ) {
	        @ids = ();
    	        ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others, $id, $cdsLength) = ("", "", "", 0, 0, 0, "", "", "", "", "", 0);
	        ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others) = split "\t", $line;
	        $fastaname =~ s/\s+.+$//;

	        if ( $attributes ne "" && $CDSfeatures eq "CDS" && $source ne "nonsense_mediated_decay" && $source ne "transposable_element" && $source ne "pseudogene" && $source ne "manual_curation_te" && $source ne "tRNA_pseudogene" && $source ne "noncoding_exon" && $source ne "mobile_element" && $source ne "misc_feature" && $source ne "gap" && $source ne "chromosome" && $source !~ /^.+RNA.*/ ) {

                if ( $firstPosition == 0 ) { $firstPosition = 1;  }
	        $cdsLength = ($secondPosition - $firstPosition) + 1;

	        #if ( $source eq "protein_coding" && $CDSfeatures eq "CDS" && $cdsLength >= 15 ) {
                    if ( $cdsLength >= 15 ) {

	                @ids = split ";", $attributes;
	                $id = $ids[0];
	                $id =~ s/^\s+//; $id =~ s/\s+$//;

	                if ( $id =~ /^ID=PAC:/ )						    { $id =~ s/^ID=PAC:/PACid:/; $id =~ s/\.+.+$//; }
	                if ( $id =~ /^ID=CDS:/ && $id =~ /\.\d+$/ )				    { $id =~ s/^ID=CDS://; $id =~ s/\.\d+$//; }
	                if ( $id =~ /^ID=CDS:/ && $id =~ /_\d+$/ )				    { $id =~ s/^ID=CDS://; $id =~ s/_\d+$//; }
	                if ( $id =~ /^ID="CDS:/ && $id =~ /_\d+"$/ )			  	    { $id =~ s/^ID="CDS://; $id =~ s/_\d+"$//; }
		        $id =~ s/^ID=CDS://;
		        $id =~ s/^ID=CDS:PAC://;
	                if ( $id =~ /^ID=/ && $id =~ /-RA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-RA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /-PA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-PA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /-TA.*$/ )				    	    { $id =~ s/^ID=//; $id =~ s/-TA.*$//; }
	                if ( $id =~ /^ID=/ && $id =~ /(:exon:|:cds:)\d+$/ )		    	    { $id =~ s/^ID=//; $id =~ s/:exon:\d+$//; $id =~ s/:cds:\d+$//; }
	                if ( $id =~ /^CDS\s+/ )						    	    { $id =~ s/^CDS\s+//; $id =~ s/\s+$//; }
	                if ( $id =~ /^ID=exon_/ )						    { $id =~ s/^ID=exon_//; $id =~ s/-\d+$//; }
	                $id =~ s/^ID=cds_of_//;
		        $id =~ s/^ID=gene://; $id =~ s/^ID=mRNA://;
	                if ( $id =~ /^ID=jCVI\|cds_/ || $id =~ /^ID=gb\|cds_/ || $id =~ /^ID=cds_/ ){ $id =~ s/^ID=jCVI\|cds_//; $id =~ s/^ID=gb\|cds_//; $id =~ s/^ID=cds_//; $id =~ s/-\d+$//; }
	                if ( $id =~ /^(ID=Trichomonas_vaginalis_Carlton\|cds_)/ )		    { $id =~ s/$1//; $id =~ s/-\d+$//; }
	            	if ( $id =~ /_exon\d+$/ )						    { $id =~ s/_exon\d+$//; $id =~ s/^ID=//; }
	            	if ( $id =~ /_exon_\d+$/ )						    { $id =~ s/_exon_\d+$//; $id =~ s/^ID=//; }
	            	if ( $id =~ /^Parent=/ && $id =~ /\.mRNA\d+$/ )			    	    { $id =~ s/^Parent=//; $id =~ s/\.mRNA\d+$//; }
	            	if ( $id =~ /^Parent=/ && $id =~ /-RA.*$/ )				    { $id =~ s/^Parent=//; $id =~ s/-RA$//; }
		    	if ( $id =~ /^Parent=/ && $id =~ /-PA.*$/ )				    { $id =~ s/^Parent=//; $id =~ s/-PA$//; }
		    	if ( $id =~ /^Parent=/ && $id =~ /-TA.*$/ )				    { $id =~ s/^Parent=//; $id =~ s/-TA$//; }
	            	if ( $id =~ /^GenePrediction\s+/ && $id =~ /-TA.*$/ )			    { $id =~ s/^GenePrediction\s+//; $id =~ s/-TA$//; }

	            	$id =~ s/^Name=//; $id =~ s/^name "//; $id =~ s/^gene_id "//; $id =~ s/"+//gis; $id =~ s/^mRNA\s+//gis; $id =~ s/^mRNA://;
	            	$id =~ s/^Name=CDS;Parent=//; $id =~ s/^Parent=//; $id =~ s/^ID=//; $id =~ s/^\|cds_//; $id =~ s/\s+.+$//;
		    	$id =~ s/(-RA|-PA):(exon|cds):\d+$//; $id =~ s/-R[A-Z]$//; $id =~ s/-P[A-Z]$//; $id =~ s/-T[A-Z]$//;

		    	$id =~ s/\.t\d+$//; $id =~ s/\.t\d+\.cds.*$//; $id =~ s/\.t\d+\.CDS.*$//;
		    	$id =~ s/\.mRNA\d+$//; $id =~ s/-mRNA-\d+$//;
		    	$id =~ s/_cds_\d+$//; $id =~ s/\.CDS\d+$//; $id =~ s/:cds$//;

	           	#print "IN SUB ENSEMBL = $db -> $source = $fastaname = $id -> $CDSfeatures --> $plusORminus:$firstPosition:$secondPosition ---> $cdsLength\n";
	            	$hashByCoordinates{$fastaname}{$plusORminus}{$CDSfeatures}{$id}{"$firstPosition:$secondPosition"} = "$firstPosition:$secondPosition";
		    	$hashByCoordinatesCheck{$fastaname}++;
		    	$cds{$id}++;
		    }
                }
	    }
	}
	close GTFF;
    }

    $cds = keys %cds;
    print "Total number of protein-coding genes parsed: $cds\n";

    return \%hashByCoordinates, \%hashByCoordinatesCheck;

    %hashByCoordinates = ();
    ($fastaname, $source, $CDSfeatures, $firstPosition, $secondPosition, $score, $plusORminus, $frame, $attributes, $others, $id, $cdsLength, $cds, $db, $gtffFile) = ("", "", "", 0, 0, 0, "", "", "", "", "", 0, 0, "", "");
}



############################################################
#                                                          #
#   DIVIDE THE FULL GENOME SEQUENCE BY FASTA FORMAT &      #
#   CALCUTE THE TOTAL NUCLEOTIDE LENGTH OF EACH SEQUENCE   #
#                                                          # 
############################################################


sub cutFastaSequences {
    # INPUT: $ref_hashByFasta = &cutFastaSequences("$pwdSequence/$sequenceFile");

    my $fastaFile = ""; my $header = ""; my $nucleotides = ""; my $fasta = ""; my $seqs = 0;
    my %hashByFasta = ();

    $fastaFile = shift;
    open SEQUENCE, "$fastaFile" || die "Cannot open sequence file\n";
    print "$fastaFile\n";

    while ( <SEQUENCE> ) {
	$fasta = $_;
        chomp $fasta;

	if ( $fasta =~ /^>/ ) {
	    $header = $fasta;
	    $header =~ s/[\s|\t]+.+$//;
	    $header =~ s/^>+[a-z]{3}\.//;
	    $header =~ s/^>+\.//;
	    $header =~ s/^>+//;
	    $header =~ s/^\s+|\s+$//gis;
	    $header =~ s/^\t+|\t+$//gis;
	    #print "$header\n";
            $seqs++;
	}

	elsif ( $fasta =~ /^[GgCcTtAaNnXx\s]/ ) {
	    $nucleotides = uc($fasta);
	    $nucleotides =~ s/\n//gis;
            $nucleotides =~ s/\s+//gis;
            $nucleotides =~ s/\t+//gis;
	    #print "$header ->\n$nucleotides\n";
	    $hashByFasta{$header} .= $nucleotides;
        }

	elsif ( $fasta !~ /^$/ && $fasta !~ /^>/ && $fasta !~ /^[GgCcTtAaNnXx\s]/) {
            print "***$fasta***\n";
	    die " ERROR: File '$fastaFile' doesn't look like a FASTA format and/or a DNA sequence!\n Please correct that and try again.";
	    exit;
	}

    }
    close SEQUENCE;

    return \%hashByFasta;
    %hashByFasta = (); $fastaFile = ""; $header = ""; $nucleotides = ""; $fasta = "";
}




############################################################################################################
#                                                                                                          #
#   IDENTIFYING THE EXON AND INTRON REFERENCE SETS AS WELL AS THEIR INTERGENIC REGIONS                     #
#                                                                                                          #
#   - getExonsIntrons        ---> obtain exon & intron coordinates and their corresponding sequence,       # 
#                                 their lengths, nucleotide frequency and intron type of length modulo 3   #
#   - getIntergenicRegions   ---> obtain the intergenic coordinates and their corresponding sequences,     #
#                                 their lengths and nucleotide frequency                                   #
#                                                                                                          # 
############################################################################################################


sub getExonsIntrons {
    # INPUT: &getExonsIntrons($minexonsize, $minintronsize, $species, "$id|location=$seqname", $gffout, $strand, \$$FastaByChromosomes{$seqname}, \@coordinatesCDS);

    my ($minexon, $minintron, $organism, $filter_header, $output, $plusORminus, $fastaSequence, $complementExon, $complementIntron, $cdsCoordinates) = ("", "", "", "", "", "", "", "", "", "");
    my ($exonStart, $exonEnd, $geneLength, $totalExonSize, $totalIntronSize) = (0, 0, 0, 0, 0);
    my ($exonSeq, $intronSeq, $multipleIntron, $pair) = ("", "", "", "");
    my ($exonSize, $intronSize, $intronStart, $intronEnd, $portionGene, $exonAT, $exonGC, $intronAT, $intronGC) = (0, 0, 0, 0, 0, 0, 0, 0, 0);
    my ($Aexon, $Texon, $Gexon, $Cexon, $Nexon, $Aintron, $Tintron, $Gintron, $Cintron, $Nintron, $ATexon, $GCexon, $ATintron, $GCintron) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    my %hashByExon = (); my %hashByIntron = (); my %hashByATGCexons = ();  my %hashByATGCintrons = ();
    my $ref_hashByATGCexon = {}; my $ref_hashByATGCintron = {};
    my @positions= (); my @starts = (); my @ends = ();
    my $i = 0; my $j = 0; my $y = 1; my $z = 1; my $x = 0;

    @positions = @{$_[7]};
    ($minexon, $minintron, $organism, $filter_header, $output, $plusORminus, $fastaSequence) = @_;
    #print "$organism, $filter_header, $output, $plusORminus, $$fastaSequence === @positions\n";

    open INTRONS, ">>$output.introns.fasta" || die "Cannot open introns.out file\n";
    open EXONS, ">>$output.exons.fasta" || die "Cannot open exons.out file\n";
    open FILTER_EXONS, ">>$output.filtered.exons.fasta" || die "Cannot open exons.out file\n";
    open FILTER_INTRONS, ">>$output.filtered.introns.fasta" || die "Cannot open introns.out file\n";


    foreach $pair ( @positions ) { ($starts[$i++], $ends[$j++]) = split ":", $pair; }
    @positions = (); 

    @starts = sort { $a <=> $b} @starts;
    @ends = sort { $a <=> $b} @ends;
    $exonStart = $starts[0];
    $exonEnd = $ends[-1];
    $cdsCoordinates = "$exonStart:$exonEnd";
    $geneLength = ($exonEnd - $exonStart) + 1;
    #print "$cdsCoordinates = $plusORminus* | *$geneLength* = *@starts* -> *@ends*\n";

    $hashByATGCexons{A} = 0;    $hashByATGCexons{T} = 0;   $hashByATGCexons{G} = 0;      $hashByATGCexons{C} = 0;    $hashByATGCexons{N} = 0;
    $hashByATGCintrons{A} = 0;  $hashByATGCintrons{T} = 0; $hashByATGCintrons{G} = 0;    $hashByATGCintrons{C} = 0;  $hashByATGCintrons{N} = 0;

    for ( $x = 0; $x <= $#starts; $x++ ) {
        $ref_hashByATGCexon = {}; $ref_hashByATGCintron = {};
	($exonSeq, $intronSeq, $multipleIntron, $complementExon, $complementIntron) = ("", "", "", "", "");
	($exonSize, $intronSize, $intronStart, $intronEnd, $portionGene, $exonAT, $exonGC, $intronAT, $intronGC) = (0, 0, 0, 0, 0, 0, 0, 0, 0);
	($Aexon, $Texon, $Gexon, $Cexon, $Nexon, $Aintron, $Tintron, $Gintron, $Cintron, $Nintron, $ATexon, $GCexon, $ATintron, $GCintron) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	$exonSeq = substr($$fastaSequence, ($starts[$x] - 1), (($ends[$x] - $starts[$x]) + 1));

        if ( $plusORminus eq "-" ) {
	    $complementExon = reverse($exonSeq);
	    $complementExon =~ tr/ATGCatgc/TACGtacg/;
            $exonSeq = "";
            $exonSeq = $complementExon;
        }
	($ref_hashByATGCexon, $exonSize) = &freqATGC(\$exonSeq);
	$portionGene = sprintf("%.2f",(($exonSize*100)/$geneLength));
	#print "exonSeq = $exonSeq || exonSize = $exonSize\n";

	if ( $exonSize >= $minexon ) {
	    $Aexon = sprintf("%.2f",((${$ref_hashByATGCexon}{A} * 100) / $exonSize)); $hashByATGCexons{A}+= ${$ref_hashByATGCexon}{A};
            $Texon = sprintf("%.2f",((${$ref_hashByATGCexon}{T} * 100) / $exonSize)); $hashByATGCexons{T}+= ${$ref_hashByATGCexon}{T};
	    $Gexon = sprintf("%.2f",((${$ref_hashByATGCexon}{G} * 100) / $exonSize)); $hashByATGCexons{G}+= ${$ref_hashByATGCexon}{G};
	    $Cexon = sprintf("%.2f",((${$ref_hashByATGCexon}{C} * 100) / $exonSize)); $hashByATGCexons{C}+= ${$ref_hashByATGCexon}{C};
	    $Nexon = sprintf("%.2f",((${$ref_hashByATGCexon}{N} * 100) / $exonSize)); $hashByATGCexons{N}+= ${$ref_hashByATGCexon}{N};
	    $exonAT = ${$ref_hashByATGCexon}{A} + ${$ref_hashByATGCexon}{T};
	    $exonGC = ${$ref_hashByATGCexon}{G} + ${$ref_hashByATGCexon}{C};
	    $ATexon = sprintf("%.2f",(($exonAT * 100) / $exonSize));
	    $GCexon = sprintf("%.2f",(($exonGC * 100) / $exonSize));
	    $totalExonSize += $exonSize;
	    $z = $y;

	    print EXONS ">$filter_header:$starts[$x]-$ends[$x]:$plusORminus:CDS_exon_$y |organism=$organism|cds=$geneLength(nts)|exon_length=$exonSize\_nts($portionGene%)|A=${$ref_hashByATGCexon}{A}($Aexon%):T=${$ref_hashByATGCexon}{T}($Texon%):G=${$ref_hashByATGCexon}{G}($Gexon%):C=${$ref_hashByATGCexon}{C}($Cexon%):AT=$exonAT($ATexon%):GC=$exonGC($GCexon%):N=${$ref_hashByATGCexon}{N}($Nexon%)\n$exonSeq\n";
            #print ">$filter_header:$starts[$x]-$ends[$x]:$plusORminus:CDS_exon_$y |organism=$organism|cds=$geneLength(nts)|exon_length=$exonSize\_nts($portionGene%)|A=${$ref_hashByATGCexon}{A}($Aexon%):T=${$ref_hashByATGCexon}{T}($Texon%):G=${$ref_hashByATGCexon}{G}($Gexon%):C=${$ref_hashByATGCexon}{C}($Cexon%):AT=$exonAT($ATexon%):GC=$exonGC($GCexon%):N=${$ref_hashByATGCexon}{N}($Nexon%)\n$exonSeq\n";

	    $hashByExon{$y++} = "$starts[$x]:$ends[$x]:$exonSize($portionGene%):${$ref_hashByATGCexon}{A}($Aexon%):${$ref_hashByATGCexon}{T}($Texon%):${$ref_hashByATGCexon}{G}($Gexon%):${$ref_hashByATGCexon}{C}($Cexon%):$exonAT($ATexon%):$exonGC($GCexon%):${$ref_hashByATGCexon}{N}($Nexon%):$geneLength";
	}

	elsif ( $exonSize < $minexon ) { 
	    print FILTER_EXONS ">$filter_header:$starts[$x]-$ends[$x]:$plusORminus:CDS_exon_$y |organism=$organism|cds=$geneLength(nts)|exon_length=$exonSize($portionGene%)\n$exonSeq\n";
	    $z = $y;
	    $y++;
	}


        unless ( $x == $#ends || ($starts[$x+1]-$ends[$x]) <= 1 ) {
	    $portionGene = 0;
	    $intronSeq = substr($$fastaSequence, $ends[$x], (($starts[$x+1] - $ends[$x]) - 1));
	    #$intronSize = length($intronSeq);

            if ( $plusORminus eq "-" ) {
	        $complementIntron = reverse($intronSeq);
	        $complementIntron =~ tr/ATGCatgc/TACGtacg/;
                $intronSeq = "";
                $intronSeq = $complementIntron;
            }
	    ($ref_hashByATGCintron, $intronSize) = &freqATGC(\$intronSeq);
	    $multipleIntron = &getIntronType($intronSize);
	    $intronStart = $ends[$x]+1;
	    $intronEnd = $starts[$x+1]-1;
	    $portionGene = sprintf("%.2f",(($intronSize*100)/$geneLength));
	    #print "intronSeq = $intronSeq || intronSize = $intronSize\n";

	    if ( $intronSize >= $minintron ) {
		$Aintron = sprintf("%.2f",((${$ref_hashByATGCintron}{A} * 100) / $intronSize)); $hashByATGCintrons{A}+= ${$ref_hashByATGCintron}{A};
		$Tintron = sprintf("%.2f",((${$ref_hashByATGCintron}{T} * 100) / $intronSize)); $hashByATGCintrons{T}+= ${$ref_hashByATGCintron}{T};
		$Gintron = sprintf("%.2f",((${$ref_hashByATGCintron}{G} * 100) / $intronSize)); $hashByATGCintrons{G}+= ${$ref_hashByATGCintron}{G};
		$Cintron = sprintf("%.2f",((${$ref_hashByATGCintron}{C} * 100) / $intronSize)); $hashByATGCintrons{C}+= ${$ref_hashByATGCintron}{C};
		$Nintron = sprintf("%.2f",((${$ref_hashByATGCintron}{N} * 100) / $intronSize)); $hashByATGCintrons{N}+= ${$ref_hashByATGCintron}{N};
		$intronAT = ${$ref_hashByATGCintron}{A} + ${$ref_hashByATGCintron}{T};
		$intronGC = ${$ref_hashByATGCintron}{G} + ${$ref_hashByATGCintron}{C};
		$ATintron = sprintf("%.2f",(($intronAT * 100) / $intronSize));
	        $GCintron = sprintf("%.2f",(($intronGC * 100) / $intronSize));
		$totalIntronSize += $intronSize;

		print INTRONS ">$filter_header:$intronStart-$intronEnd:$plusORminus:CDS_intron_$z:$multipleIntron |organism=$organism|cds=$geneLength(nts)|intron_length=$intronSize\_nts($portionGene%)|A=${$ref_hashByATGCintron}{A}($Aintron%):T=${$ref_hashByATGCintron}{T}($Tintron%):${$ref_hashByATGCintron}{G}($Gintron%):C=${$ref_hashByATGCintron}{C}($Cintron%):AT=$intronAT($ATintron%):GC=$intronGC($GCintron%):N=${$ref_hashByATGCintron}{N}($Nintron%)\n$intronSeq\n";
		#print ">$filter_header:$intronStart-$intronEnd:$plusORminus:CDS_intron_$z:$multipleIntron |organism=$organism|cds=$geneLength(nts)|intron_length=$intronSize\_nts($portionGene%)|A=${$ref_hashByATGCintron}{A}($Aintron%):T=${$ref_hashByATGCintron}{T}($Tintron%):${$ref_hashByATGCintron}{G}($Gintron%):C=${$ref_hashByATGCintron}{C}($Cintron%):AT=$intronAT($ATintron%):GC=$intronGC($GCintron%):N=${$ref_hashByATGCintron}{N}($Nintron%)\n$intronSeq\n";

		$hashByIntron{$z} = "$intronStart:$intronEnd:$intronSize($portionGene%):${$ref_hashByATGCintron}{A}($Aintron%):${$ref_hashByATGCintron}{T}($Tintron%):${$ref_hashByATGCintron}{G}($Gintron%):${$ref_hashByATGCintron}{C}($Cintron%):$intronAT($ATintron%):$intronGC($GCintron%):${$ref_hashByATGCintron}{N}($Nintron%):$geneLength:$multipleIntron";
	    }

	    elsif ( $intronSize < $minintron ) {
		print FILTER_INTRONS ">$filter_header:$intronStart-$intronEnd:$plusORminus:CDS_intron_$z:$multipleIntron |organism=$organism|cds=$geneLength(nts)|intron_length=$intronSize($portionGene%)\n$intronSeq\n";
		#$z++;
	    }
	}
    }

    return ($totalExonSize, $totalIntronSize, $cdsCoordinates, \%hashByExon, \%hashByIntron, \%hashByATGCexons, \%hashByATGCintrons, $geneLength);

    %hashByExon = (); %hashByIntron = (); %hashByATGCexons = (); %hashByATGCintrons = ();
    $ref_hashByATGCexon = {}; $ref_hashByATGCintron = {};
    @starts = (); @ends = ();
}


sub getIntergenicRegions {
    # INPUT: ($ref_hashByIntergenics, $ref_hashByATGCintergenics) = &getIntergenicRegions($species, "$id|location=$seqname", $gffout, $strand, $lengthSequence, \$$FastaByChromosomes{$seqname}, \@intergenics);

    my ($organism, $filter_header, $output, $plusORminus, $sizeSequence, $fastaSequence, $getFirstPosition, $getLastPosition, $lost_right, $lost_left) = ("", "", "", "", "", "", "", "", "", "");
    my ($intergenicSeq, $getSequence, $pair, ) = ("", "", "");
    my ($intergenicSize, $intergenicStart, $intergenicEnd, $portionSequence, $intergenicAT, $intergenicGC, $lastSize) = (0, 0, 0, 0, 0, 0, 0);
    my ($Aintergenic, $Tintergenic, $Gintergenic, $Cintergenic, $Nintergenic, $ATintergenic, $GCintergenic) = (0, 0, 0, 0, 0, 0, 0);
    my %hashByIntergenic = (); my %hashByATGCintergenics = (); my $ref_hashByATGCintergenic = {};
    my @positions = (); my @starts = (); my @ends = ();
    my $i = 0; my $j = 0; my $y = 1; my $z = 1; my $x = 0;

    @positions = @{$_[6]};
    ($organism, $filter_header, $output, $plusORminus, $sizeSequence, $fastaSequence) = @_;
    $hashByATGCintergenics{A} = 0; $hashByATGCintergenics{T} = 0; $hashByATGCintergenics{G} = 0; 
    $hashByATGCintergenics{C} = 0; $hashByATGCintergenics{N} = 0;

    open INTERGENICS, ">>$output.intergenics.fasta" || die "Cannot open intergenics.out file\n";

    foreach $pair ( @positions ) { ($starts[$i++], $ends[$j++]) = split ":", $pair; }
    @starts = sort { $a <=> $b} @starts;
    @ends = sort { $a <=> $b} @ends;
    #print "*$plusORminus* | *$sizeSequence* = *@starts* -> *@ends*\n";

    ($getFirstPosition, $lost_right) = split ":", $positions[0];
    ($lost_left, $getLastPosition) = split ":", $positions[-1];
    $lastSize = $sizeSequence - $getLastPosition;
    #print "1 - $getFirstPosition ==== $getLastPosition - $sizeSequence\n";
    @positions = (); 

    if ( $getFirstPosition > 1 ) {
        $ref_hashByATGCintergenic = {};
	($intergenicSeq, $getSequence) = ("", "");
	($intergenicSize, $intergenicStart, $intergenicEnd, $portionSequence, $intergenicAT, $intergenicGC) = (0, 0, 0, 0, 0, 0);
	($Aintergenic, $Tintergenic, $Gintergenic, $Cintergenic, $Nintergenic, $ATintergenic, $GCintergenic) = (0, 0, 0, 0, 0, 0, 0);

        if ( $plusORminus eq "+" ) {
	    $intergenicSeq = substr($$fastaSequence, 0, ($getFirstPosition - 1));
	}
	if ( $plusORminus eq "-" ) {
	    $getSequence = substr($$fastaSequence, 0, ($getFirstPosition - 1));
	    $intergenicSeq = reverse($getSequence);
	    $intergenicSeq =~ tr/ATGCatgc/TACGtacg/;
	}
	($ref_hashByATGCintergenic, $intergenicSize) = &freqATGC(\$intergenicSeq);
	#print "FIRST intergenicSeq = $intergenicSeq || FIRST intergenicSize = $intergenicSize\n";

	if ( $intergenicSize > 0 ) {
	    $Aintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{A} * 100) / $intergenicSize)); $hashByATGCintergenics{A} += ${$ref_hashByATGCintergenic}{A};
	    $Tintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{T} * 100) / $intergenicSize)); $hashByATGCintergenics{T} += ${$ref_hashByATGCintergenic}{T};
	    $Gintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{G} * 100) / $intergenicSize)); $hashByATGCintergenics{G} += ${$ref_hashByATGCintergenic}{G};
	    $Cintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{C} * 100) / $intergenicSize)); $hashByATGCintergenics{C} += ${$ref_hashByATGCintergenic}{C};
	    $Nintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{N} * 100) / $intergenicSize)); $hashByATGCintergenics{N} += ${$ref_hashByATGCintergenic}{N};
	    $intergenicAT = ${$ref_hashByATGCintergenic}{A} + ${$ref_hashByATGCintergenic}{T};
	    $intergenicGC = ${$ref_hashByATGCintergenic}{G} + ${$ref_hashByATGCintergenic}{C};
	    $ATintergenic = sprintf("%.2f",(($intergenicAT * 100) / $intergenicSize));
	    $GCintergenic = sprintf("%.2f",(($intergenicGC * 100) / $intergenicSize));
	    $intergenicStart = 1;
	    $intergenicEnd = $getFirstPosition - 1;
	    $portionSequence = sprintf("%.2f",(($intergenicSize*100)/$sizeSequence));

            print INTERGENICS ">$filter_header:$plusORminus:$intergenicStart-$intergenicEnd:intergenic_$y($intergenicSize nts) |organism=$species|length=$intergenicSize($portionSequence%)|A=${$ref_hashByATGCintergenic}{A}($Aintergenic%):T=${$ref_hashByATGCintergenic}{T}($Tintergenic%):G=${$ref_hashByATGCintergenic}{G}($Gintergenic%):C=${$ref_hashByATGCintergenic}{C}($Cintergenic%):AT=$intergenicAT($ATintergenic%):GC=$intergenicGC($GCintergenic%):N=${$ref_hashByATGCintergenic}{N}($Nintergenic%)\n$intergenicSeq\n";
	    #print "FIRST INTERGENIC -> $intergenicSeq=$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize\n";

	    $hashByIntergenic{$y++} = "$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize";
        }
    }

    for ( $x = 0; $x <= $#starts; $x++ ) {
        $ref_hashByATGCintergenic = {};
	($intergenicSeq, $getSequence) = ("", "");
	($intergenicSize, $intergenicStart, $intergenicEnd, $portionSequence, $intergenicAT, $intergenicGC) = (0, 0, 0, 0, 0, 0);
	($Aintergenic, $Tintergenic, $Gintergenic, $Cintergenic, $Nintergenic, $ATintergenic, $GCintergenic) = (0, 0, 0, 0, 0, 0, 0);

	unless ( $x == $#ends || ($starts[$x+1]-$ends[$x]) <= 1 ) {
	    if ( $plusORminus eq "+" ) {
	        $intergenicSeq = substr($$fastaSequence, $ends[$x], (($starts[$x+1] - $ends[$x]) - 1));
	    }
	    if ( $plusORminus eq "-" ) {
	        $getSequence = substr($$fastaSequence, $ends[$x], (($starts[$x+1] - $ends[$x]) - 1));
	        $intergenicSeq = reverse($getSequence);
	        $intergenicSeq =~ tr/ATGCatgc/TACGtacg/;
	    }
	    ($ref_hashByATGCintergenic, $intergenicSize) = &freqATGC(\$intergenicSeq);
	    #print "MIDDLE intergenicSeq = $intergenicSeq || MIDDLE intergenicSize = $intergenicSize\n";

	    if ( $intergenicSize > 0 ) {
	        $Aintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{A} * 100) / $intergenicSize)); $hashByATGCintergenics{A} += ${$ref_hashByATGCintergenic}{A};
		$Tintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{T} * 100) / $intergenicSize)); $hashByATGCintergenics{T} += ${$ref_hashByATGCintergenic}{T};
		$Gintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{G} * 100) / $intergenicSize)); $hashByATGCintergenics{G} += ${$ref_hashByATGCintergenic}{G};
		$Cintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{C} * 100) / $intergenicSize)); $hashByATGCintergenics{C} += ${$ref_hashByATGCintergenic}{C};
		$Nintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{N} * 100) / $intergenicSize)); $hashByATGCintergenics{N} += ${$ref_hashByATGCintergenic}{N};
	        $intergenicAT = ${$ref_hashByATGCintergenic}{A} + ${$ref_hashByATGCintergenic}{T};
	        $intergenicGC = ${$ref_hashByATGCintergenic}{G} + ${$ref_hashByATGCintergenic}{C};
	        $ATintergenic = sprintf("%.2f",(($intergenicAT * 100) / $intergenicSize));
	        $GCintergenic = sprintf("%.2f",(($intergenicGC * 100) / $intergenicSize));
	        $intergenicStart = $ends[$x]+1;
	        $intergenicEnd = $starts[$x+1]-1;
	        $portionSequence = sprintf("%.2f",(($intergenicSize*100)/$sizeSequence));

                print INTERGENICS ">$filter_header:$plusORminus:$intergenicStart-$intergenicEnd:intergenic_$y($intergenicSize nts) |organism=$species|length=$intergenicSize($portionSequence%)|A=${$ref_hashByATGCintergenic}{A}($Aintergenic%):T=${$ref_hashByATGCintergenic}{T}($Tintergenic%):G=${$ref_hashByATGCintergenic}{G}($Gintergenic%):C=${$ref_hashByATGCintergenic}{C}($Cintergenic%):AT=$intergenicAT($ATintergenic%):GC=$intergenicGC($GCintergenic%):N=${$ref_hashByATGCintergenic}{N}($Nintergenic%)\n$intergenicSeq\n";

	        $hashByIntergenic{$y++} = "$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize";
	        #print "MIDDLE INTERGENIC -> $intergenicSeq=$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize\n";
            }
	}
    }

    if ( $getLastPosition < $sizeSequence ) {
        $ref_hashByATGCintergenic = {};
    	($intergenicSeq, $getSequence) = ("", "");
	($intergenicSize, $intergenicStart, $intergenicEnd, $portionSequence, $intergenicAT, $intergenicGC) = (0, 0, 0, 0, 0, 0);
	($Aintergenic, $Tintergenic, $Gintergenic, $Cintergenic, $Nintergenic, $ATintergenic, $GCintergenic) = (0, 0, 0, 0, 0, 0, 0);

        if ( $plusORminus eq "+" ) {
	    $intergenicSeq = substr($$fastaSequence, $getLastPosition, $lastSize);
	}
	if ( $plusORminus eq "-" ) {
	    $getSequence = substr($$fastaSequence, $getLastPosition, $lastSize);
	    $intergenicSeq = reverse($getSequence);
	    $intergenicSeq =~ tr/ATGCatgc/TACGtacg/;
	}
	($ref_hashByATGCintergenic, $intergenicSize) = &freqATGC(\$intergenicSeq);
	#print "LAST intergenicSeq = $intergenicSeq || LAST intergenicSize = $intergenicSize\n";

	if ( $intergenicSize > 0 ) {
	    $Aintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{A} * 100) / $intergenicSize)); $hashByATGCintergenics{A} += ${$ref_hashByATGCintergenic}{A};
	    $Tintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{T} * 100) / $intergenicSize)); $hashByATGCintergenics{T} += ${$ref_hashByATGCintergenic}{T};
	    $Gintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{G} * 100) / $intergenicSize)); $hashByATGCintergenics{G} += ${$ref_hashByATGCintergenic}{G};
	    $Cintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{C} * 100) / $intergenicSize)); $hashByATGCintergenics{C} += ${$ref_hashByATGCintergenic}{C};
	    $Nintergenic = sprintf("%.2f",((${$ref_hashByATGCintergenic}{N} * 100) / $intergenicSize)); $hashByATGCintergenics{N} += ${$ref_hashByATGCintergenic}{N};
	    $intergenicAT = ${$ref_hashByATGCintergenic}{A} + ${$ref_hashByATGCintergenic}{T};
	    $intergenicGC = ${$ref_hashByATGCintergenic}{G} + ${$ref_hashByATGCintergenic}{C};
	    $ATintergenic = sprintf("%.2f",(($intergenicAT * 100) / $intergenicSize));
	    $GCintergenic = sprintf("%.2f",(($intergenicGC * 100) / $intergenicSize));
	    $intergenicStart = $getLastPosition + 1;
	    $intergenicEnd = $sizeSequence;
	    $portionSequence = sprintf("%.2f",(($intergenicSize*100)/$sizeSequence));

            print INTERGENICS ">$filter_header:$plusORminus:$intergenicStart-$intergenicEnd:intergenic_$y($intergenicSize nts) |organism=$species|length=$intergenicSize($portionSequence%)|A=${$ref_hashByATGCintergenic}{A}($Aintergenic%):T=${$ref_hashByATGCintergenic}{T}($Tintergenic%):G=${$ref_hashByATGCintergenic}{G}($Gintergenic%):C=${$ref_hashByATGCintergenic}{C}($Cintergenic%):AT=$intergenicAT($ATintergenic%):GC=$intergenicGC($GCintergenic%):N=${$ref_hashByATGCintergenic}{N}($Nintergenic%)\n$intergenicSeq\n";

	    $hashByIntergenic{$y++} = "$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize";
	    #print "LAST INTERGENIC -> $intergenicSeq=$intergenicStart:$intergenicEnd:$intergenicSize($portionSequence%):${$ref_hashByATGCintergenic}{A}($Aintergenic%):${$ref_hashByATGCintergenic}{T}($Tintergenic%):${$ref_hashByATGCintergenic}{G}($Gintergenic%):${$ref_hashByATGCintergenic}{C}($Cintergenic%):$intergenicAT($ATintergenic%):$intergenicGC($GCintergenic%):${$ref_hashByATGCintergenic}{N}($Nintergenic%):$intergenicSize\n";
        }
    }

    return (\%hashByIntergenic, \%hashByATGCintergenics);
    %hashByIntergenic = (); %hashByATGCintergenics = (); $ref_hashByATGCintergenic = {}; $intergenicSeq = ""; $getSequence = "";
    @starts = (); @ends = ();
}




######################################################################################################
#                                                                                                    #
#   CALCULATING SEVERAL STATISTIC INDIVIDUAL DESCRIPTORS                                             #
#                                                                                                    #
#   - getStatisticsBySpecies ---> calculate the frequency of nucleotides in a sequence                      #
#   - getIntronType   ---> estimate the intron type of length modulo 3                               #
#   - getFreqATGCtype ---> calculate the fraction of nucleotide types from the total sequence size   #
#                                                                                                    # 
######################################################################################################


sub getStatisticsBySpecies {
    my $type = ""; my $sps = "";     my $sequence = ""; my $size_sp = 0; my $cds_sp = 0;
    my $dir = "";  my $process = ""; my $outfile_gff = ""; my $output = "";   my $definition = "";

    my %coordinates   = ();  my %cds_chr        = ();  my %cds_str          = ();  my %number_sp     = ();  my %number_chr     = ();  my %number_str = ();
    my %length_sp     = ();  my %length_chr     = ();  my %length_str       = ();  my %mean_sp       = ();  my %mean_chr       = ();  my %mean_str = ();
    my %atgc_sp       = ();  my %atgc_chr       = ();  my %atgc_str         = ();  my %wavg_sp       = ();  my %wavg_chr       = ();  my %wavg_str = (); 
    my %statistics_sp = ();  my %statistics_chr = ();  my %statistics_str   = ();  my %intronsCDS_sp = ();  my %intronsCDS_chr = ();  my %intronsCDS_str = ();
    my %nts_sp        = ();  my %size_chr       = ();  my %ncdnasize_chr    = ();  my %ncprotein_chr = ();  my %density_sp     = ();
    my %all_outfile   = ();  my %order_exons    = ();  my %orderexons_stats = ();  my %order_introns = ();  my %orderintrons_stats = ();

    my @frequency = (); my @definition = ();

    $type = shift; $sps = shift; $sequence = shift; $dir = shift; $process = shift; $outfile_gff = shift; $size_sp = shift; $cds_sp = shift;
    %coordinates   = %{$_[0]};   %cds_chr        = %{$_[1]};   %cds_str        = %{$_[2]};
    %number_sp     = %{$_[3]};   %number_chr     = %{$_[4]};   %number_str     = %{$_[5]};
    %length_sp     = %{$_[6]};   %length_chr     = %{$_[7]};   %length_str     = %{$_[8]};
    %atgc_sp       = %{$_[9]};   %atgc_chr       = %{$_[10]};  %atgc_str       = %{$_[11]}; 
    %statistics_sp = %{$_[12]};  %statistics_chr = %{$_[13]};  %statistics_str = %{$_[14]};
    %nts_sp	   = %{$_[15]};  %size_chr       = %{$_[16]};  %all_outfile    = %{$_[17]};

    if ( $type eq "exons" || $type eq "introns" ) {
        %mean_sp = %{$_[18]}; %mean_chr = %{$_[19]}; %mean_str = %{$_[20]}; %wavg_sp = %{$_[21]}; %wavg_chr = %{$_[22]}; %wavg_str = %{$_[23]}; %density_sp = %{$_[24]};
    }
    if ( $type eq "introns" ) {
        %order_exons = %{$_[25]}; %orderexons_stats = %{$_[26]}; %order_introns = %{$_[27]}; %orderintrons_stats = %{$_[28]};
	%intronsCDS_sp = %{$_[29]}; %intronsCDS_chr = %{$_[30]}; %intronsCDS_str = %{$_[31]}; @definition = @{$_[32]}; $definition = pop;
    }
    if ( $type eq "intergenics" ) {
        %ncprotein_chr = %{$_[18]}; %ncdnasize_chr = %{$_[19]};
    }

    $output = $outfile_gff;
    $outfile_gff .= ".coordinates.txt";
    #print "$output.exons.content.txt\n";

    if ( $process =~ /[NnFf]/ ) {
        open STATSEXONS_CHROMOSOME, ">>$dir/total.exons.chromosomes.table.txt" || die "Cannot open $dir/total.exons.chromosomes.table.txt\n";
        open STATSEXONS_STRANDS, ">>$dir/total.exons.strands.table.txt" || die "Cannot open $dir/total.exons.strands.table.txt\n";
        open STATSEXONS_GENOME, ">>$dir/total.exons.table.txt" || die "Cannot open $dir/total.exons.table.txt\n";

        open STATSINTRONS_CHROMOSOME, ">>$dir/total.introns.chromosomes.table.txt" || die "Cannot open $dir/total.introns.chromosomes.table.txt\n";
        open STATSINTRONS_STRANDS, ">>$dir/total.introns.strands.table.txt" || die "Cannot open $dir/total.introns.strands.table.txt\n";
        open STATSINTRONS_GENOME, ">>$dir/total.introns.table.txt" || die "Cannot open $dir/total.introns.table.txt\n";

        open STATSINTERGENICS_CHROMOSOME, ">>$dir/total.intergenics.chromosomes.table.txt" || die "Cannot open $dir/total.intergenics.chromosomes.table.txt\n";
        open STATSINTERGENICS_STRANDS, ">>$dir/total.intergenics.strands.table.txt" || die "Cannot open $dir/total.intergenics.strands.table.txt\n";
        open STATSINTERGENICS_GENOME, ">>$dir/total.intergenics.table.txt" || die "Cannot open $dir/total.intergenics.table.txt\n";

        open EXONSRANGE_GENOME, ">>$dir/total.exons_ranges.table.txt" || die "Cannot open $dir/total.exons_ranges.table.txt\n";
        open INTRONSRANGE_GENOME, ">>$dir/total.introns_ranges.table.txt" || die "Cannot open $dir/total.introns_ranges.table.txt\n";
        open INTERGENICSRANGE_GENOME, ">>$dir/total.intergenics_ranges.table.txt" || die "Cannot open $dir/total.intergenics_ranges.table.txt\n";

        open EXONSATGC_GENOME, ">>$dir/total.exons_ranges_atgc.table.txt" || die "Cannot open $dir/total.exons_ranges_atgc.table.txt\n";
        open INTRONSATGC_GENOME, ">>$dir/total.introns_ranges_atgc.table.txt" || die "Cannot open $dir/total.introns_ranges_atgc.table.txt\n";
        open INTERGENICSATGC_GENOME, ">>$dir/total.intergenics_ranges_atgc.table.txt" || die "Cannot open $dir/total.intergenics_ranges_atgc.table.txt\n";

	open SIZESEXONS_FILTER, ">>$output.quartiles.sizes.exons.distribution.genome.txt" || die "Cannot open $output.exons.content.txt\n";
	open SIZESINTRONS_FILTER, ">>$output.quartiles.sizes.introns.distribution.genome.txt" || die "Cannot open $output.introns.content.txt\n";
	open SIZESINTERGENICS_FILTER, ">>$output.quartiles.sizes.intergenics.distribution.genome.txt" || die "Cannot open $output.intergenics.content.txt\n";

	open SIZESEXONS_ALL, ">>$output.all.sizes.exons.distribution.genome.txt" || die "Cannot open $output.exons.content.txt\n";
	open SIZESINTRONS_ALL, ">>$output.all.sizes.introns.distribution.genome.txt" || die "Cannot open $output.introns.content.txt\n";
	open SIZESINTERGENICS_ALL, ">>$output.all.sizes.intergenics.distribution.genome.txt" || die "Cannot open $output.intergenics.content.txt\n";

	open COORDINATES_EXONS, ">>$output.exons.content.txt" || die "Cannot open $output.exons.content.txt\n";
	open COORDINATES_INTRONS, ">>$output.introns.content.txt" || die "Cannot open $output.introns.content.txt\n";
	open COORDINATES_INTERGENICS, ">>$output.intergenics.content.txt" || die "Cannot open $output.intergenics.content.txt\n";
	
	open DEFINITION, ">>$output.order.intron-exon-intron.distribution.txt" || die "Cannot open $output.order.intron-exon-intron.distribution.txt\n";
	open ORDER_EXONS, ">>$output.order.exons.distribution.txt" || die "Cannot open $output.order.exons.distribution.txt\n";
	open ORDER_INTRONS, ">>$output.order.introns.distribution.txt" || die "Cannot open $output.order.introns.distribution.txt\n";

	open DEFINITION_FILTER, ">>$output.quartiles.order.intron-exon-intron.distribution.txt" || die "Cannot open $output.order.intron-exon-intron.distribution.txt\n";
	open ORDEREXONS_FILTER, ">>$output.quartiles.order.exons.distribution.txt" || die "Cannot open $output.order.exons.distribution.txt\n";
	open ORDERINTRONS_FILTER, ">>$output.quartiles.order.introns.distribution.txt" || die "Cannot open $output.order.introns.distribution.txt\n";
    }
    if ( $process =~ /[YyTt]/ ) {
        open STATSEXONS_CHROMOSOME, ">>$dir/total.$sps.exons.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.exons.chromosomes.table.txt\n";
        open STATSEXONS_STRANDS, ">>$dir/total.$sps.exons.strands.table.txt" || die "Cannot open $dir/total.$sps.exons.strands.table.txt\n";
        open STATSEXONS_GENOME, ">>$dir/total.$sps.exons.table.txt" || die "Cannot open $dir/total.$sps.exons.table.txt\n";

        open STATSINTRONS_CHROMOSOME, ">>$dir/total.$sps.introns.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.introns.chromosomes.table.txt\n";
        open STATSINTRONS_STRANDS, ">>$dir/total.$sps.introns.strands.table.txt" || die "Cannot open $dir/total.$sps.introns.strands.table.txt\n";
        open STATSINTRONS_GENOME, ">>$dir/total.$sps.introns.table.txt" || die "Cannot open $dir/total.$sps.introns.table.txt\n";

        open STATSINTERGENICS_CHROMOSOME, ">>$dir/total.$sps.intergenics.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.intergenics.chromosomes.table.txt\n";
        open STATSINTERGENICS_STRANDS, ">>$dir/total.$sps.intergenics.strands.table.txt" || die "Cannot open $dir/total.$sps.intergenics.strands.table.txt\n";
        open STATSINTERGENICS_GENOME, ">>$dir/total.$sps.intergenics.table.txt" || die "Cannot open $dir/total.$sps.intergenics.table.txt\n";

        open EXONSRANGE_GENOME, ">>$dir/total.$sps.exons_ranges.table.txt" || die "Cannot open $dir/total.$sps.exons_ranges.table.txt\n";
        open INTRONSRANGE_GENOME, ">>$dir/total.$sps.introns_ranges.table.txt" || die "Cannot open $dir/total.$sps.introns_ranges.table.txt\n";
        open INTERGENICSRANGE_GENOME, ">>$dir/total.$sps.intergenics_ranges.table.txt" || die "Cannot open $dir/total.$sps.intergenics_ranges.table.txt\n";

        open EXONSATGC_GENOME, ">>$dir/total.$sps.exons_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.exons_ranges_atgc.table.txt\n";
        open INTRONSATGC_GENOME, ">>$dir/total.$sps.introns_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.introns_ranges_atgc.table.txt\n";
        open INTERGENICSATGC_GENOME, ">>$dir/total.$sps.intergenics_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.intergenics_ranges_atgc.table.txt\n";

	open SIZESEXONS_FILTER, ">>$output.quartiles.sizes.exons.distribution.genome.txt" || die "Cannot open $output.exons.content.txt\n";
	open SIZESINTRONS_FILTER, ">>$output.quartiles.sizes.introns.distribution.genome.txt" || die "Cannot open $output.introns.content.txt\n";
	open SIZESINTERGENICS_FILTER, ">>$output.quartiles.sizes.intergenics.distribution.genome.txt" || die "Cannot open $output.intergenics.content.txt\n";

	open SIZESEXONS_ALL, ">>$output.all.sizes.exons.distribution.genome.txt" || die "Cannot open $output.exons.content.txt\n";
	open SIZESINTRONS_ALL, ">>$output.all.sizes.introns.distribution.genome.txt" || die "Cannot open $output.introns.content.txt\n";
	open SIZESINTERGENICS_ALL, ">>$output.all.sizes.intergenics.distribution.genome.txt" || die "Cannot open $output.intergenics.content.txt\n";

	open COORDINATES_EXONS, ">>$output.exons.content.txt" || die "Cannot open $output.exons.content.txt\n";
	open COORDINATES_INTRONS, ">>$output.introns.content.txt" || die "Cannot open $output.introns.content.txt\n";
	open COORDINATES_INTERGENICS, ">>$output.intergenics.content.txt" || die "Cannot open $output.intergenics.content.txt\n";

	open DEFINITION, ">>$output.order.intron-exon-intron.distribution.txt" || die "Cannot open $output.order.intron-exon-intron.distribution.txt\n";
	open ORDER_EXONS, ">>$output.order.exons.distribution.txt" || die "Cannot open $output.order.exons.distribution.txt\n";
	open ORDER_INTRONS, ">>$output.order.introns.distribution.txt" || die "Cannot open $output.order.introns.distribution.txt\n";

	open DEFINITION_FILTER, ">>$output.quartiles.order.intron-exon-intron.distribution.txt" || die "Cannot open $output.order.intron-exon-intron.distribution.txt\n";
	open ORDEREXONS_FILTER, ">>$output.quartiles.order.exons.distribution.txt" || die "Cannot open $output.order.exons.distribution.txt\n";
	open ORDERINTRONS_FILTER, ">>$output.quartiles.order.introns.distribution.txt" || die "Cannot open $output.order.introns.distribution.txt\n";
    }


    my %Ranges = ();                 my %RangesATGC = ();           my %SizesFrequency = ();       my %DensityFrequency = (); 
    my $ref_hashDistribution = {};   my $ref_hashRangesATGC = {};   my $ref_hashSizesCount = {};   my $ref_hashDensityCount = {};

    my $r10 = 0;        my $r50 = 0;         my $r100 = 0;        my $r250 = 0;         my $r500 = 0;         my $r1000 = 0;
    my $r5000 = 0;      my $r10000 = 0;      my $r50000 = 0;      my $r100000 = 0;      my $r500000 = 0;      my $r1000000 = 0;
    my $r5000000 = 0;   my $r10000000 = 0;   my $r50000000 = 0;   my $r100000000 = 0;   my $r500000000 = 0;   my $r1000000000 = 0; 

    my @y_10 = ();       my @y_50 = ();        my @y_100 = ();       my @y_250 = ();        my @y_500 = ();        my @y_1000 = ();
    my @y_5000 = ();     my @y_10000 = ();     my @y_50000 = ();     my @y_100000 = ();     my @y_500000 = ();     my @y_1000000 = ();
    my @y_5000000 = ();  my @y_10000000 = ();  my @y_50000000 = ();  my @y_100000000 = ();  my @y_500000000 = ();  my @y_1000000000 = ();

    my @OverlapsExons       = ();  my @LengthsOverlapsExons       = ();  my @coordinates         = ();  my @density           = ();  
    my @OverlapsIntrons     = ();  my @LengthsOverlapsIntrons     = ();  my @IntersectsIntrons   = ();  my @atgcs             = (); 
    my @OverlapsIntergenics = ();  my @LengthsOverlapsIntergenics = ();  my @exonpos_line        = ();  my @intronerrors      = ();
    my @OverlapsCDS         = ();  my @LengthsIntersectsIntrons   = ();  my @intronpositions     = ();  my @means             = ();
    my @densitypos_filter   = ();  my @density_filter             = ();  my @densitypos_all      = ();  my @density_all       = ();
    my @definition_all      = ();  my @definition_filter          = ();  my @pairs_quartiles     = ();  my @pairs_all         = (); 
    my @exonerrors_filter   = ();  my @intronave_filter           = ();  my @intronerrors_filter = ();  my @intronpos_filter  = (); 
    my @exonave_filter      = ();  my @exonpos_filter             = ();  my @exon_counts         = ();  my @intron_counts     = ();
    my @exonaverages        = ();  my @exonpositions              = ();  my @exonerrors          = ();  my @intronaverages    = ();
    my @avg_density         = ();  my @avg_densityfilter          = ();  my @weighted_density    = ();  my @ratio_density     = ();
    my @statistics          = ();  my @intronpos_line             = ();  my @high_density        = ();  my @lower_ratio       = ();
    my @intronsizes         = ();  my @exonsizes                  = ();  my @x_ranges            = ();  my @filter            = ();
    my @LengthsOverlapsCDS  = ();  my $ref_arrayOverlapsIntrons   = ();  my $ref_arrayOverlaps   = ();
    my $ref_arrayLength     = ();  my $ref_arrayLengthIntrons     = ();

    my $TranscriptStructure = "";       my $IntronDensityExons = "";  my $coordinates = "";           my $sizes = "";
    my $avg_densityfilter = "";         my $seqname = "";             my $strand = "";                my $keys = "";
    my $weighted_density = "";          my $chart_density = "";       my $x_all = "";                 my $x_filter = "";
    my $chart3d_definition = "";        my $avg_density = "";         my $dataSet = "";               my $filter = "";
    my $line_introns = "";              my $line_exons = "";          my $x_introns = "";             my $x_exons = "";
    my $chart_errors = "";              my $chart_order = "";         my $exon_errors = "";           my $intron_errors = "";
    my $ratio_density = "";             my $par = "";                 my $sum = "";                   my $num = "";

    my $avg_filter = 0;                 my $wavg_filter = 0;          my $totalSizeIntron = 0;        my $percentage = 0;
    my $median = 0;                     my $null = 0;                 my $average = 0;                my $waverage = 0;
    my $variance = 0;                   my $stddev = 0;               my $var_filter = 0;             my $sd_filter = 0;
    my $variance_wavg = 0;              my $stddev_wavg = 0;          my $varwavg_filter = 0;         my $sdwavg_filter = 0;
    my $errorbarBegin = 0;              my $errorbarEnd = 0;          my $errorbarBeginFilter = 0;    my $errorbarEndFilter = 0; 
    my $errorbarBegin_wavg = 0;         my $errorbarEnd_wavg = 0;     my $exonNumber = 0;             my $intronNumber = 0;
    my $errorbarBeginFilter_wavg = 0;   my $logGenome = 0;            my $totalSize = 0;              my $totalCDS = 0;
    my $errorbarEndFilter_wavg = 0;     my $logChromosome = 0;        my $logSizeChromosome = 0;      my $logMeanChromosome = 0;
    my $GeneDensityChromosome = 0;      my $MeanSizeChromosome = 0;   my $ChromosomeFrequency = 0;    my $chart_ranges = 0;
    my $AverageSizeChromosome = 0;      my $LengthMbsChromosome = 0;  my $par1 = 0;                   my $par2 = 0;
    my $GenomeFrequencyTotalCDS = 0;    my $logSizeGenome = 0;        my $logAverageGenome = 0;       my $logMeanGenome = 0;
    my $TotalNTSchromosome = 0;         my $ATfreqChromosome = 0;     my $GCfreqChromosome = 0;       my $NfreqChromosome = 0; 
    my $AfreqChromosome = 0;            my $TfreqChromosome = 0;      my $GfreqChromosome = 0;        my $CfreqChromosome = 0;
    my $logAverageChromosome = 0;       my $logSizeStrand = 0;        my $logAverageStrand = 0;       my $logMeanStrand = 0;
    my $OverlapsIntergenics = 0;        my $OverlapsExons = 0;        my $IntersectsIntrons = 0;      my $exon_shortintron = 0;
    my $GenomeFrequencyMbsCDS = 0;      my $GenomeFrequencyCDS = 0;   my $WeightedAverageGenome = 0;  my $GenomeFrequencyTotal = 0;
    my $GeneDensityGenome = 0;          my $AverageSizeGenome = 0;    my $LengthMbsGenome = 0;        my $GenomeFrequencyMbs = 0;
    my $GenomeFrequency = 0;            my $MeanSizeGenome = 0;       my $intronsCDSgenome = 0;       my $FilterAverageGenome = 0;
    my $TotalNTSgenome = 0;             my $ATfreqGenome = 0;         my $GCfreqGenome = 0;           my $NfreqGenome = 0;
    my $AfreqGenome = 0;                my $TfreqGenome = 0;          my $GfreqGenome = 0;            my $CfreqGenome = 0;
    my $intronsCDSchromosome = 0;       my $variance_wden = 0;        my $stddev_wden = 0;            my $varwden_filter = 0;
    my $ChromosomeFrequencyCDS = 0;     my $wavgden_filter = 0;       my $daverage = 0;               my $sdwden_filter = 0;
    my $errorbarBeginFilter_wden = 0;   my $def_exons = 0;            my $def_total = 0;              my $deffilter_total = 0;
    my $errorbarEndFilter_wden = 0;     my $StrandFrequency = 0;      my $MeanSizeStrand = 0;         my $intronsCDSstrand = 0;
    my $errorbarBegin_wden = 0;         my $errorbarEnd_wden = 0;     my $summatory = 0;              my $Dupper_fence = 0;
    my $exonSizeByPosAve_filter = 0;    my $SequenceSizeMbs = 0;      my $exon_filterintron = 0;      my $Dlower_fence = 0;
    my $exonSizeByPositionAve = 0;      my $defexons_filter = 0;      my $x_densityfilter = 0;        my $x_density = 0;
    my $intronSizeByPositionAve = 0;    my $density_filter = 0;       my $scalar = 0;                 my $Dmedian = 0;  
    my $GenomeSizeMbs = 0;              my $GeneDensityStrand = 0;    my $AverageSizeStrand = 0;      my $LengthMbsStrand = 0; 
    my $AfreqStrand = 0;                my $TfreqStrand = 0;          my $GfreqStrand = 0;            my $CfreqStrand = 0;
    my $TotalNTSstrand = 0;             my $ATfreqStrand = 0;         my $GCfreqStrand = 0;           my $NfreqStrand = 0;
    my $intronSizeByPosAve_filter = 0;
    my $LengthMbsOverlapChromosomeCDS = 0; 
    my $LengthMbsOverlapChromosome = 0;

    my $Q1 = 0;  my $Q3 = 0;  my $Q_top = 0;  my $Q_down = 0;  my $IQR = 0; my $lower_fence = 0;  my $upper_fence = 0;  my $x = 0;
    my $DQ1 = 0; my $DQ3 = 0; my $DQ_top = 0; my $DQ_down = 0; my $DIQR = 0; my $freqGC = 0; my $freqAT = 0; my $freqNN = 0;


    for ( $x = 0; $x <= 17 ; $x++ ) {
	push @y_exons, "0";
        push @y_introns, "0";
        push @y_intergenics, "0";
    }

    $GenomeSizeMbs = sprintf("%.2f",($size_sp / 1000000));

    if ( defined @{ $nts_sp{$sps}{$type} } ) {
        @atgcs = sort {$a cmp $b} @{ $nts_sp{$sps}{$type} };
        ($ref_hashDistribution, $ref_hashRangesATGC, $summatory) = &getTypeSizeDistribution(\@atgcs);
        %Ranges = %{$ref_hashDistribution};
        %RangesATGC = %{$ref_hashRangesATGC};
	#print "$type = @atgcs\n";
    }
    else { $freqGC = 0; $freqAT = 0; $freqNN = 0; }


    if ( $type eq "exons" ) {
        @y_exons = (); $sum_exons = 0; $sum_exons = "($summatory)";
        print EXONSRANGE_GENOME "$sps ($summatory)\t$GenomeSizeMbs\t";

        for $keys ( sort {$a<=>$b} keys %Ranges ) {
	    $percentage = 0; $freqGC = 0; $freqAT = 0; $freqNN = 0;

	    if ( $Ranges{$keys} > 0 ) { $percentage = sprintf("%.2f",(($Ranges{$keys} * 100) / $summatory)); }
	    print EXONSRANGE_GENOME "$percentage\t";
	    #print "$keys = $percentage\n";
	    push @y_exons, "$percentage";

	    if ( $RangesATGC{$keys}{GC} > 0 ) { $freqGC = sprintf("%.2f",(($RangesATGC{$keys}{GC} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{AT} > 0 ) { $freqAT = sprintf("%.2f",(($RangesATGC{$keys}{AT} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{N} > 0 )  { $freqNN = sprintf("%.2f",(($RangesATGC{$keys}{N} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    print EXONSATGC_GENOME "$keys\t$Ranges{$keys}\t$freqGC\t$freqAT\t$freqNN\n";
	}
	print EXONSRANGE_GENOME "\n";


        $GenomeFrequencyTotal = 0;
        for $seqname ( sort keys %{ $size_chr{$sps} } ) {
            @statistics = (); @means = ();
	    $logChromosome = 0; $logSizeChromosome = 0; $logAverageChromosome = 0; $logMeanChromosome = 0;
            $coordinates = ""; $OverlapsExons = 0; $SequenceSizeMbs = 0; $GeneDensityChromosome = 0;
	    $AverageSizeChromosome = 0; $LengthMbsChromosome = 0; $LengthMbsOverlapChromosome = 0;
	    $ATfreqChromosome = 0; $GCfreqChromosome = 0; $NfreqChromosome = 0; $AfreqChromosome = 0; 
            $TfreqChromosome = 0; $GfreqChromosome = 0; $CfreqChromosome = 0; $TotalNTSchromosome = 0; 
	    $MeanSizeChromosome = 0; $LengthMbsOverlapChromosome = 0; $ChromosomeFrequency = 0;
	    $Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0;
            $median = 0; $null = ""; $average = 0; $waverage = 0;
	    $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $variance_wavg = 0;
            $stddev_wavg = 0; $varwavg_filter = 0; $sdwavg_filter = 0; $avg_filter = 0; $wavg_filter = 0;
	    $errorbarBegin = 0; $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0; 
            $errorbarBegin_wavg = 0; $errorbarEnd_wavg = 0; $errorbarBeginFilter_wavg = 0; $errorbarEndFilter_wavg = 0;

	    if ( $cds_chr{$sps}{$seqname} > 0 && $number_chr{$sps}{$seqname}{$type} > 0 ) {
	        $SequenceSizeMbs = sprintf("%.2f",($size_chr{$sps}{$seqname} / 1000000));
	        $GeneDensityChromosome = sprintf("%.2f",($number_chr{$sps}{$seqname}{$type} / $cds_chr{$sps}{$seqname}));

### average length of all exons in the genome:
	        $AverageSizeChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / $number_chr{$sps}{$seqname}{$type}));
                $LengthMbsChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / 1000000));

### weighted average length of exons in CDS:
	        $MeanSizeChromosome = sprintf("%.2f",($mean_chr{$sps}{$seqname}{$type} / $cds_chr{$sps}{$seqname}));		

		($TotalNTSchromosome, $AfreqChromosome, $TfreqChromosome, $GfreqChromosome, $CfreqChromosome, $NfreqChromosome, $ATfreqChromosome, $GCfreqChromosome) = getFreqATGCtype("chromosomes", $sps, $seqname, $type, \%atgc_chr);

	        @OverlapsExons = (); @LengthsOverlapsExons = (); @coordinates = (); $totalSize = 0; $ref_arrayOverlaps = (); $ref_arrayLength = ();
	        @coordinates = sort {$a cmp $b} @{ $coordinates{$sps}{$seqname}{exons} };
	        ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getOverlaps("$size_chr{$sps}{$seqname}", \@coordinates);
	        @OverlapsExons = @$ref_arrayOverlaps;
	        @LengthsOverlapsExons = @$ref_arrayLength;
	        $LengthMbsOverlapChromosome = sprintf("%.2f",($totalSize / 1000000));
                $ChromosomeFrequency = sprintf("%.2f",(($totalSize * 100) / $size_chr{$sps}{$seqname}));
	        $GenomeFrequencyTotal+= $totalSize;
		$coordinates = scalar @coordinates;
		$OverlapsExons = scalar @OverlapsExons;
	        #print "*$seqname*\t$#coordinates = @coordinates\n$#OverlapsExons = @OverlapsExons\n-$totalSize-$LengthMbsOverlapChromosome\n";
		#print ">$sps:$seqname|$size_chr{$sps}{$seqname}_nts|exons:$#coordinates:$#OverlapsExons|$totalSize\_nts:$LengthMbsOverlapChromosome\_Mbs:$ChromosomeFrequency%\n@OverlapsExons\n";

		print COORDINATES_EXONS ">$sps:$seqname|$size_chr{$sps}{$seqname}_nts|exons:$coordinates:$OverlapsExons|$totalSize\_nts:$LengthMbsOverlapChromosome\_Mbs:$ChromosomeFrequency%\n@OverlapsExons\n";

	        @statistics = sort {$a <=> $b} @{ $statistics_chr{$sps}{$seqname}{exons} };
		($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);
	        ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($AverageSizeChromosome, $lower_fence, $upper_fence, \@statistics);
		$errorbarBegin = sprintf("%.2f",($AverageSizeChromosome - $stddev)); $errorbarEnd = sprintf("%.2f",($AverageSizeChromosome + $stddev));
		$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));

		@means = sort {$a <=> $b} @{ $wavg_chr{$sps}{$seqname}{$type} };
		($variance_wavg, $stddev_wavg, $varwavg_filter, $sdwavg_filter, $wavg_filter, $null, $waverage) = &calculateStdDev($MeanSizeChromosome, $lower_fence, $upper_fence, \@means);
		$errorbarBegin_wavg = sprintf("%.2f",($MeanSizeChromosome - $stddev_wavg)); $errorbarEnd_wavg = sprintf("%.2f",($MeanSizeChromosome + $stddev_wavg));
		$errorbarBeginFilter_wavg = sprintf("%.2f",($wavg_filter - $sdwavg_filter)); $errorbarEndFilter_wavg = sprintf("%.2f",($wavg_filter + $sdwavg_filter));

		if ( $SequenceSizeMbs > 0 ) { $logChromosome = sprintf("%.4f", (log10($SequenceSizeMbs))); } else { $logChromosome = 0; }
		if ( $LengthMbsOverlapChromosome > 0 ) { $logSizeChromosome = sprintf("%.4f", (log10($LengthMbsOverlapChromosome))); } else { $logSizeChromosome = 0; }
	    	if ( $AverageSizeChromosome > 0 ) { $logAverageChromosome = sprintf("%.4f", (log10($AverageSizeChromosome))); } else { $logAverageChromosome = 0; }
	    	if ( $MeanSizeChromosome > 0 ) { $logMeanChromosome = sprintf("%.4f", (log10($MeanSizeChromosome))); } else { $logMeanChromosome = 0; }

	        if ( $sequence eq "genome" ) {
                    print STATSEXONS_CHROMOSOME "$sps\t$seqname\t$SequenceSizeMbs\t$cds_chr{$sps}{$seqname}\t$GeneDensityChromosome\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$MeanSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\t$logMeanChromosome\n";
                    #print "STATSEXONS_CHROMOSOME: $sps\t$seqname\t$SequenceSizeMbs\t$cds_chr{$sps}{$seqname}\t$GeneDensityChromosome\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$MeanSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\t$logMeanChromosome\n";
		}
	    }
	}
    }

    if ( $type eq "introns" ) {
        @y_introns = (); $sum_introns = 0; $sum_introns = "($summatory)";
        print INTRONSRANGE_GENOME "$sps ($summatory)\t$GenomeSizeMbs\t";

        for $keys ( sort {$a<=>$b} keys %Ranges ) {
	    $percentage = 0; $freqGC = 0; $freqAT = 0; $freqNN = 0;

	    if ( $Ranges{$keys} > 0 ) { $percentage = sprintf("%.2f",(($Ranges{$keys} * 100) / $summatory)); }
	    print INTRONSRANGE_GENOME "$percentage\t";
	    push @y_introns, "$percentage";

	    if ( $RangesATGC{$keys}{GC} > 0 ) { $freqGC = sprintf("%.2f",(($RangesATGC{$keys}{GC} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{AT} > 0 ) { $freqAT = sprintf("%.2f",(($RangesATGC{$keys}{AT} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{N} > 0 )  { $freqNN = sprintf("%.2f",(($RangesATGC{$keys}{N} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    print INTRONSATGC_GENOME "$keys\t$Ranges{$keys}\t$freqGC\t$freqAT\t$freqNN\n";
	}
	print INTRONSRANGE_GENOME "\n";


        $GenomeFrequencyTotal = 0;
        for $seqname ( sort keys %{ $size_chr{$sps} } ) {
            @statistics = (); @means = ();
	    $logChromosome = 0; $logSizeChromosome = 0; $logAverageChromosome = 0; $logMeanChromosome = 0;
            $coordinates = ""; $IntersectsIntrons = 0; $SequenceSizeMbs = 0; $GeneDensityChromosome = 0;
            $AverageSizeChromosome = 0; $LengthMbsChromosome = 0; $LengthMbsOverlapChromosome = 0;
	    $TotalNTSchromosome = 0; $ATfreqChromosome = 0; $GCfreqChromosome = 0; $NfreqChromosome = 0;
            $AfreqChromosome = 0; $TfreqChromosome = 0; $GfreqChromosome = 0; $CfreqChromosome = 0;
	    $MeanSizeChromosome = 0; $intronsCDSchromosome = 0; $LengthMbsOverlapChromosome = 0; $ChromosomeFrequency = 0;
	    $Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0;
            $median = 0; $null = ""; $average = 0; $waverage = 0;
	    $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $variance_wavg = 0;
            $stddev_wavg = 0; $varwavg_filter = 0; $sdwavg_filter = 0; $avg_filter = 0; $wavg_filter = 0;
	    $errorbarBegin = 0; $errorbarEnd = 0; $errorbarBegin_wavg = 0; $errorbarEnd_wavg = 0;
            $errorbarBeginFilter = 0; $errorbarEndFilter = 0; $errorbarBeginFilter_wavg = 0; $errorbarEndFilter_wavg = 0;

	    if ( $cds_chr{$sps}{$seqname} > 0 && $number_chr{$sps}{$seqname}{$type} > 0 ) {
	        $SequenceSizeMbs = sprintf("%.2f",($size_chr{$sps}{$seqname} / 1000000));
	        $GeneDensityChromosome = sprintf("%.2f",($number_chr{$sps}{$seqname}{$type} / $intronsCDS_chr{$sps}{$seqname}));

### average length of all introns in the genome:
	        $AverageSizeChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / $number_chr{$sps}{$seqname}{$type}));
                $LengthMbsChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / 1000000));

### weighted average length of introns in CDS bearing introns:
	        $MeanSizeChromosome = sprintf("%.2f",($mean_chr{$sps}{$seqname}{$type} / $intronsCDS_chr{$sps}{$seqname}));		
                $intronsCDSchromosome = sprintf("%.2f",(($intronsCDS_chr{$sps}{$seqname} * 100) / $cds_chr{$sps}{$seqname}));

		($TotalNTSchromosome, $AfreqChromosome, $TfreqChromosome, $GfreqChromosome, $CfreqChromosome, $NfreqChromosome, $ATfreqChromosome, $GCfreqChromosome) = getFreqATGCtype("chromosomes", $sps, $seqname, $type, \%atgc_chr);

	        @OverlapsExons = (); @LengthsOverlapsExons = (); @coordinates = (); $totalSize = 0; $ref_arrayOverlaps = (); $ref_arrayLength = ();
	        @coordinates = sort {$a cmp $b} @{ $coordinates{$sps}{$seqname}{exons} };
	        ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getOverlaps("$size_chr{$sps}{$seqname}", \@coordinates);
	        @OverlapsExons = @$ref_arrayOverlaps;
	        @LengthsOverlapsExons = @$ref_arrayLength;
	        #print "*$seqname*\t$#coordinates = @coordinates\n$#OverlapsExons = @OverlapsExons\n-$totalSize-\n";
	        #print "*$seqname* -> overlaped exons = $#coordinates -> $#OverlapsExons\n";

	        @OverlapsIntrons = (); @LengthsOverlapsIntrons = (); @coordinates = (); $totalSize = 0; $ref_arrayOverlaps = (); $ref_arrayLength = ();
	        @coordinates = sort {$a cmp $b} @{ $coordinates{$sps}{$seqname}{introns} };
	        ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getOverlaps("$size_chr{$sps}{$seqname}", \@coordinates);
	        @OverlapsIntrons = @$ref_arrayOverlaps;
	        @LengthsOverlapsIntrons = @$ref_arrayLength;
	        #print "*$seqname*\t$#coordinates = @coordinates\n$#OverlapsIntrons = @OverlapsIntrons\n-$totalSize-\n";
		#print "*$seqname* -> overlaped introns = $#coordinates -> $#OverlapsIntrons\n";

### -------------------- ###
		@IntersectsIntrons = (); @LengthsIntersectsIntrons = (); $totalSizeIntron = 0; $ref_arrayOverlapsIntrons = (); $ref_arrayLengthIntrons = ();
	        ($totalSizeIntron, $ref_arrayOverlapsIntrons, $ref_arrayLengthIntrons) = &getIntersections("$size_chr{$sps}{$seqname}", \@OverlapsIntrons, \@OverlapsExons);
		@IntersectsIntrons = @$ref_arrayOverlapsIntrons;
		@LengthsIntersectsIntrons = @$ref_arrayLengthIntrons;
	        $LengthMbsOverlapChromosome = sprintf("%.2f",($totalSizeIntron / 1000000));
                $ChromosomeFrequency = sprintf("%.2f",(($totalSizeIntron * 100) / $size_chr{$sps}{$seqname}));
	        $GenomeFrequencyTotal+= $totalSizeIntron;
		$coordinates = scalar @coordinates;
		$IntersectsIntrons = scalar @IntersectsIntrons;
		print COORDINATES_INTRONS ">$sps:$seqname|$size_chr{$sps}{$seqname}_nts|introns:$coordinates:$IntersectsIntrons|$totalSizeIntron\_nts:$LengthMbsOverlapChromosome\_Mbs:$ChromosomeFrequency%\n@IntersectsIntrons\n";
### -------------------- ###

	        @statistics = sort {$a <=> $b} @{ $statistics_chr{$sps}{$seqname}{$type} };
		($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);
		($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($AverageSizeChromosome, $lower_fence, $upper_fence, \@statistics);
		$errorbarBegin = sprintf("%.2f",($AverageSizeChromosome - $stddev)); $errorbarEnd = sprintf("%.2f",($AverageSizeChromosome + $stddev));
		$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));

		@means = sort {$a <=> $b} @{ $wavg_chr{$sps}{$seqname}{$type} };
		($variance_wavg, $stddev_wavg, $varwavg_filter, $sdwavg_filter, $wavg_filter, $null, $waverage) = &calculateStdDev($MeanSizeChromosome, $lower_fence, $upper_fence, \@means);
		$errorbarBegin_wavg = sprintf("%.2f",($MeanSizeChromosome - $stddev_wavg)); $errorbarEnd_wavg = sprintf("%.2f",($MeanSizeChromosome + $stddev_wavg));
		$errorbarBeginFilter_wavg = sprintf("%.2f",($wavg_filter - $sdwavg_filter)); $errorbarEndFilter_wavg = sprintf("%.2f",($wavg_filter + $sdwavg_filter));

		if ( $SequenceSizeMbs > 0 ) { $logChromosome = sprintf("%.4f", (log10($SequenceSizeMbs))); } else { $logChromosome = 0; }
		if ( $LengthMbsOverlapChromosome > 0 ) { $logSizeChromosome = sprintf("%.4f", (log10($LengthMbsOverlapChromosome))); } else { $logSizeChromosome = 0; }
	    	if ( $AverageSizeChromosome > 0 ) { $logAverageChromosome = sprintf("%.4f", (log10($AverageSizeChromosome))); } else { $logAverageChromosome = 0; }
	    	if ( $MeanSizeChromosome > 0 ) { $logMeanChromosome = sprintf("%.4f", (log10($MeanSizeChromosome))); } else { $logMeanChromosome = 0; }

	        if ( $sequence eq "genome" ) {
                    print STATSINTRONS_CHROMOSOME "$sps\t$seqname\t$SequenceSizeMbs\t$cds_chr{$sps}{$seqname}\t$intronsCDSchromosome\t$GeneDensityChromosome\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$MeanSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\t$logMeanChromosome\n";
                    #print "STATSINTRONS_CHROMOSOME: $sps\t$seqname\t$SequenceSizeMbs\t$cds_chr{$sps}{$seqname}\t$intronsCDSchromosome\t$GeneDensityChromosome\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$MeanSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\t$logMeanChromosome\n";
		}
	    }
	}
    }

    if ( $type eq "intergenics" ) {
        @y_intergenics = (); $sum_intergenics = 0; $sum_intergenics = "($summatory)";
        print INTERGENICSRANGE_GENOME "$sps ($summatory)\t$GenomeSizeMbs\t";

        for $keys ( sort {$a<=>$b} keys %Ranges ) {
	    $percentage = 0; $freqGC = 0; $freqAT = 0; $freqNN = 0;

	    if ( $Ranges{$keys} > 0 ) { $percentage = sprintf("%.2f",(($Ranges{$keys} * 100) / $summatory)); }
	    print INTERGENICSRANGE_GENOME "$percentage\t";
	    push @y_intergenics, "$percentage";

	    if ( $RangesATGC{$keys}{GC} > 0 ) { $freqGC = sprintf("%.2f",(($RangesATGC{$keys}{GC} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{AT} > 0 ) { $freqAT = sprintf("%.2f",(($RangesATGC{$keys}{AT} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    if ( $RangesATGC{$keys}{N} > 0 )  { $freqNN = sprintf("%.2f",(($RangesATGC{$keys}{N} * 100) / ($RangesATGC{$keys}{AT} + $RangesATGC{$keys}{GC} + $RangesATGC{$keys}{N}))); }
	    print INTERGENICSATGC_GENOME "$keys\t$Ranges{$keys}\t$freqGC\t$freqAT\t$freqNN\n";
	}
	print INTERGENICSRANGE_GENOME "\n";


        $GenomeFrequencyTotal = 0;
	$GenomeFrequencyTotalCDS = 0;
        for $seqname ( sort keys %{ $size_chr{$sps} } ) {
            @statistics = ();
	    $logChromosome = 0; $logSizeChromosome = 0; $logAverageChromosome = 0; $coordinates = ""; $OverlapsIntergenics = 0;
	    $SequenceSizeMbs = 0; $AverageSizeChromosome = 0; $LengthMbsChromosome = 0; $LengthMbsOverlapChromosome = 0;
	    $TotalNTSchromosome = 0; $ATfreqChromosome = 0; $GCfreqChromosome = 0; $NfreqChromosome = 0;
            $AfreqChromosome = 0; $TfreqChromosome = 0; $GfreqChromosome = 0; $CfreqChromosome = 0;
	    $ChromosomeFrequencyCDS = 0; $LengthMbsOverlapChromosomeCDS = 0; $LengthMbsOverlapChromosome = 0; $ChromosomeFrequency = 0;
	    $Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0; $median = 0; $null = ""; $average = 0;
	    $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $avg_filter = 0;
	    $errorbarBeginFilter = 0; $errorbarEndFilter = 0; $errorbarBeginFilter_wavg = 0; $errorbarEndFilter_wavg = 0;

	    if ( $cds_chr{$sps}{$seqname} > 0 && $number_chr{$sps}{$seqname}{$type} > 0 ) {
	        $SequenceSizeMbs = sprintf("%.2f",($size_chr{$sps}{$seqname} / 1000000));
	        $AverageSizeChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / $number_chr{$sps}{$seqname}{$type}));
                $LengthMbsChromosome = sprintf("%.2f",($length_chr{$sps}{$seqname}{$type} / 1000000));

		($TotalNTSchromosome, $AfreqChromosome, $TfreqChromosome, $GfreqChromosome, $CfreqChromosome, $NfreqChromosome, $ATfreqChromosome, $GCfreqChromosome) = getFreqATGCtype("chromosomes", $sps, $seqname, $type, \%atgc_chr);

	        @OverlapsCDS = (); @LengthsOverlapsCDS = (); @coordinates = (); $totalCDS = ""; $ref_arrayOverlaps = (); $ref_arrayLength = ();
	        @coordinates = sort {$a cmp $b} @{ $coordinates{$sps}{$seqname}{intergenics} };
	        ($totalCDS, $ref_arrayOverlaps, $ref_arrayLength) = &getOverlaps("$size_chr{$sps}{$seqname}", \@coordinates);
	        @OverlapsCDS = @$ref_arrayOverlaps;
	        @LengthsOverlapsCDS = @$ref_arrayLength;
	        #print "*$seqname*\t$#coordinates = @coordinates\n$#OverlapsExons = @OverlapsExons\n-$totalSize-\n";
		$LengthMbsOverlapChromosomeCDS = sprintf("%.2f",($totalCDS / 1000000));
                $ChromosomeFrequencyCDS = sprintf("%.2f",(($totalCDS * 100) / $size_chr{$sps}{$seqname}));
	        $GenomeFrequencyTotalCDS+= $totalCDS;

	        @OverlapsIntergenics = (); @LengthsOverlapsIntergenics = (); $totalSize = 0; $ref_arrayOverlaps = (); $ref_arrayLength = ();
	        ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getNonCDSoverlaps("$size_chr{$sps}{$seqname}", \@coordinates);
	        @OverlapsIntergenics = @$ref_arrayOverlaps;
	        @LengthsOverlapsIntergenics = @$ref_arrayLength;
	        #print "*$seqname*\t$#coordinates = @coordinates\n$#OverlapsExons = @OverlapsExons\n-$totalSize-\n";
	        $LengthMbsOverlapChromosome = sprintf("%.2f",($totalSize / 1000000));
                $ChromosomeFrequency = sprintf("%.2f",(($totalSize * 100) / $size_chr{$sps}{$seqname}));
	        $GenomeFrequencyTotal+= $totalSize;
		$coordinates = scalar @coordinates;
		$OverlapsIntergenics = scalar @OverlapsIntergenics;
		print COORDINATES_INTERGENICS ">$sps:$seqname|$size_chr{$sps}{$seqname}_nts|intergenics:$OverlapsIntergenics|$totalSize\_nts:$LengthMbsOverlapChromosome\_Mbs:$ChromosomeFrequency%\n@OverlapsIntergenics\n";

	        @statistics = sort {$a <=> $b} @{ $statistics_chr{$sps}{$seqname}{$type} };
		($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);
		($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($AverageSizeChromosome, $lower_fence, $upper_fence, \@statistics);
		$errorbarBegin = sprintf("%.2f",($AverageSizeChromosome - $stddev)); $errorbarEnd = sprintf("%.2f",($AverageSizeChromosome + $stddev));
		$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));

		if ( $SequenceSizeMbs > 0 ) { $logChromosome = sprintf("%.4f", (log10($SequenceSizeMbs))); } else { $logChromosome = 0; }
		if ( $LengthMbsOverlapChromosome > 0 ) { $logSizeChromosome = sprintf("%.4f", (log10($LengthMbsOverlapChromosome))); } else { $logSizeChromosome = 0; }
	    	if ( $AverageSizeChromosome > 0 ) { $logAverageChromosome = sprintf("%.4f", (log10($AverageSizeChromosome))); } else { $logAverageChromosome = 0; }

	        if ( $sequence eq "genome" ) {
	            print STATSINTERGENICS_CHROMOSOME "$sps\t$seqname\t$SequenceSizeMbs\t$cds_sp\t$LengthMbsOverlapChromosomeCDS\t$ChromosomeFrequencyCDS\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\n";
	            #print "STATSINTERGENICS_CHROMOSOME: $sps\t$seqname\t$SequenceSizeMbs\t$cds_sp\t$LengthMbsOverlapChromosomeCDS\t$ChromosomeFrequencyCDS\t$number_chr{$sps}{$seqname}{$type}\t$LengthMbsChromosome\t$LengthMbsOverlapChromosome\t$ChromosomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeChromosome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqChromosome\t$TfreqChromosome\t$GfreqChromosome\t$CfreqChromosome\t$NfreqChromosome\t$ATfreqChromosome\t$GCfreqChromosome\t$logChromosome\t$logSizeChromosome\t$logAverageChromosome\n";
		}
	    }
	    else { print COORDINATES_INTERGENICS ">$sps:$seqname|$size_chr{$sps}{$seqname}_nts|intergenics:1|$ncdnasize_chr{$sps}{$seqname}{nts}\_nts:$ncdnasize_chr{$sps}{$seqname}{mbs}\_Mbs:100%\n1:$ncdnasize_chr{$sps}{$seqname}{nts}\n"; }
	}
	$GenomeFrequencyTotal+= $ncprotein_chr{$sps};
    }
    print "**$type -> chromosomes finished**\n";

    for $strand ( sort keys %{ $number_str{$sps}{$type} } ) {
        @statistics = (); @means = ();
	$logGenome = 0; $logSizeStrand = 0; $logAverageStrand = 0; $logMeanStrand = 0;
	$GenomeSizeMbs = 0; $GeneDensityStrand = 0; $AverageSizeStrand = 0; $LengthMbsStrand = 0;
        $StrandFrequency = 0; $MeanSizeStrand = 0; $intronsCDSstrand = 0; $TotalNTSstrand = 0;
	$NfreqStrand = 0; $AfreqStrand = 0; $TfreqStrand = 0; $GfreqStrand = 0; $CfreqStrand = 0;
	$Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0;
        $median = 0; $null = ""; $average = 0; $waverage = 0; $ATfreqStrand = 0; $GCfreqStrand = 0;
	$variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0;
        $variance_wavg = 0; $stddev_wavg = 0; $varwavg_filter = 0; $sdwavg_filter = 0;
        $avg_filter = 0; $wavg_filter = 0; $errorbarBegin = 0; $errorbarEnd = 0;
        $errorbarBegin_wavg = 0; $errorbarEnd_wavg = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;
        $errorbarBeginFilter_wavg = 0; $errorbarEndFilter_wavg = 0;

	if ( $cds_str{$sps}{$strand} > 0 && $number_str{$sps}{$type}{$strand} > 0 ) {
	    $GenomeSizeMbs = sprintf("%.2f",($size_sp / 1000000));

### standard average length of all features in a chromosome-strand = 
### total length of all features in a strand of a chromosome divided by the total number of all features in the strand of that chromosome:

            $AverageSizeStrand = sprintf("%.2f",($length_str{$sps}{$type}{$strand} / $number_str{$sps}{$type}{$strand}));	
            $LengthMbsStrand = sprintf("%.2f",($length_str{$sps}{$type}{$strand} / 1000000));
	    $StrandFrequency = sprintf("%.2f",(($length_str{$sps}{$type}{$strand} * 100) / $size_sp));

	    if ( $type eq "introns" ) {
	        $intronsCDSstrand = sprintf("%.2f",(($intronsCDS_str{$sps}{$strand} * 100) / $cds_str{$sps}{$strand}));
		$GeneDensityStrand = sprintf("%.2f",($number_str{$sps}{$type}{$strand} / $intronsCDS_str{$sps}{$strand}));

### weighted average length of introns in CDS bearing introns:
		$MeanSizeStrand = sprintf("%.2f",($mean_str{$sps}{$type}{$strand} / $intronsCDS_str{$sps}{$strand}));
	    }
	    if ( $type eq "exons" ) {
	        $GeneDensityStrand = sprintf("%.2f",($number_str{$sps}{$type}{$strand} / $cds_str{$sps}{$strand}));

### weighted average length of exons in CDS:
	        $MeanSizeStrand = sprintf("%.2f",($mean_str{$sps}{$type}{$strand} / $cds_str{$sps}{$strand})); 			
	    }

	    ($TotalNTSstrand, $AfreqStrand, $TfreqStrand, $GfreqStrand, $CfreqStrand, $NfreqStrand, $ATfreqStrand, $GCfreqStrand) = getFreqATGCtype($strand, $sps, $seqname, $type, \%atgc_str);
	    @statistics = sort {$a <=> $b} @{ $statistics_str{$sps}{$type}{$strand} };
	    ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);
	    ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($AverageSizeStrand, $lower_fence, $upper_fence, \@statistics);
	    $errorbarBegin = sprintf("%.2f",($AverageSizeStrand - $stddev)); $errorbarEnd = sprintf("%.2f",($AverageSizeStrand + $stddev));
	    $errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));

	    if ( $type eq "exons" || $type eq "introns" ) {
	        @means = sort {$a <=> $b} @{ $wavg_str{$sps}{$type}{$strand} };
	        ($variance_wavg, $stddev_wavg, $varwavg_filter, $sdwavg_filter, $wavg_filter, $null, $waverage) = &calculateStdDev($MeanSizeStrand, $lower_fence, $upper_fence, \@means);
		$errorbarBegin_wavg = sprintf("%.2f",($MeanSizeStrand - $stddev_wavg)); $errorbarEnd_wavg = sprintf("%.2f",($MeanSizeStrand + $stddev_wavg));
	    	$errorbarBeginFilter_wavg = sprintf("%.2f",($wavg_filter - $sdwavg_filter)); $errorbarEndFilter_wavg = sprintf("%.2f",($wavg_filter + $sdwavg_filter));
	    }

	    $logGenome = sprintf("%.4f", (log10($GenomeSizeMbs)));
	    if ( $LengthMbsStrand > 0 ) { $logSizeStrand = sprintf("%.4f", (log10($LengthMbsStrand))); } else { $logSizeStrand = 0; }
	    if ( $AverageSizeStrand > 0 ) { $logAverageStrand = sprintf("%.4f", (log10($AverageSizeStrand))); } else { $logAverageStrand = 0; }
	    if ( $MeanSizeStrand > 0 ) { $logMeanStrand = sprintf("%.4f", (log10($MeanSizeStrand))); } else { $logMeanStrand = 0; }

	    if ( $type eq "exons" ) {
                print STATSEXONS_STRANDS "$sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$GeneDensityStrand\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$MeanSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\t$logMeanStrand\n";
                #print "STATSEXONS_STRANDS: $sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$GeneDensityStrand\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$MeanSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\t$logMeanStrand\n";
            }
	    if ( $type eq "introns" ) {
                print STATSINTRONS_STRANDS "$sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$intronsCDSstrand\t$GeneDensityStrand\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$MeanSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\t$logMeanStrand\n";
                #print "STATSINTRONS_STRANDS: $sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$intronsCDSstrand\t$GeneDensityStrand\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$MeanSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\t$logMeanStrand\n";
	    }
	    if ( $type eq "intergenics" ) {
	        print STATSINTERGENICS_STRANDS "$sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\n";
	        #print "STATSINTERGENICS_STRANDS: $sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t$number_str{$sps}{$type}{$strand}\t$LengthMbsStrand\t$StrandFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeStrand\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqStrand\t$TfreqStrand\t$GfreqStrand\t$CfreqStrand\t$NfreqStrand\t$ATfreqStrand\t$GCfreqStrand\t$logGenome\t$logSizeStrand\t$logAverageStrand\n";
            }
	}
	if ( $type eq "introns" && !defined $number_str{$sps}{introns}{$strand} ) {
            print STATSINTRONS_STRANDS "$sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
            #print "STATSINTRONS_STRANDS: $sps\t$GenomeSizeMbs\t$strand\t$cds_str{$sps}{$strand}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
	}
    }
    print "**$type -> strands finished**\n";

    if ( $cds_sp > 0 && $number_sp{$sps}{$type} > 0 ) {
        @statistics = (); %SizesFrequency = (); %DensityFrequency = (); @means = (); @density = ();
	$logGenome = 0; $logSizeGenome = 0; $logAverageGenome = 0; $logMeanGenome = 0;
        $GenomeFrequencyMbsCDS = 0; $GenomeFrequencyCDS = 0; $WeightedAverageGenome = 0;
        $TranscriptStructure = 0; $IntronDensityExons = 0; $GenomeSizeMbs = 0; $GeneDensityGenome = 0;
        $AverageSizeGenome = 0; $LengthMbsGenome = 0; $GenomeFrequencyMbs = 0; $GenomeFrequency = 0;
        $MeanSizeGenome = 0; $intronsCDSgenome = 0; $FilterAverageGenome = 0; $TotalNTSgenome = 0;
        $ATfreqGenome = 0; $GCfreqGenome = 0; $NfreqGenome = 0; $AfreqGenome = 0; $TfreqGenome = 0;
        $GfreqGenome = 0; $CfreqGenome = 0; $Dupper_fence = 0; $Dmedian = 0; 
	$Q1 = 0; $Q3 = 0; $Q_top = 0; $Q_down = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0;
        $median = 0; $DQ1 = 0; $DQ3 = 0; $DQ_top = 0; $DQ_down = 0; $DIQR = 0; $Dlower_fence = 0; 
	$variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $variance_wavg = 0; $stddev_wavg = 0;
        $varwavg_filter = 0; $sdwavg_filter = 0; $avg_filter = 0; $wavg_filter = 0;
	$variance_wden = 0; $stddev_wden = 0; $varwden_filter = 0; $sdwden_filter = 0; $wavgden_filter = 0;
        $null = ""; $average = 0; $waverage = 0; $daverage = 0;
	$errorbarBegin = 0; $errorbarEnd = 0; $errorbarBegin_wavg = 0; $errorbarEnd_wavg = 0; $errorbarBeginFilter = 0;
        $errorbarEndFilter = 0; $errorbarBeginFilter_wavg = 0; $errorbarEndFilter_wavg = 0;
	$errorbarBegin_wden = 0; $errorbarEnd_wden = 0; $errorbarBeginFilter_wden = 0; $errorbarEndFilter_wden = 0;
        $wavgden_filter = 0; $ref_hashSizesCount = {}; $ref_hashDensityCount = {};

        $GenomeFrequencyMbs = sprintf("%.2f",($GenomeFrequencyTotal / 1000000));
        $GenomeFrequency = sprintf("%.2f",(($GenomeFrequencyTotal * 100) / $size_sp));
        $GenomeFrequencyMbsCDS = sprintf("%.2f",($GenomeFrequencyTotalCDS / 1000000));
        $GenomeFrequencyCDS = sprintf("%.2f",(($GenomeFrequencyTotalCDS * 100) / $size_sp));
        $GenomeSizeMbs = sprintf("%.2f",($size_sp / 1000000));

### average length of all features in the genome:
        $AverageSizeGenome = sprintf("%.2f",($length_sp{$sps}{$type} / $number_sp{$sps}{$type}));
        $LengthMbsGenome = sprintf("%.2f",($length_sp{$sps}{$type} / 1000000));

	($TotalNTSgenome, $AfreqGenome, $TfreqGenome, $GfreqGenome, $CfreqGenome, $NfreqGenome, $ATfreqGenome, $GCfreqGenome) = getFreqATGCtype("genome", $sps, $seqname, $type, \%atgc_sp);

        @statistics = sort {$a <=> $b} @{ $statistics_sp{$sps}{$type} };
	($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);
	($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $ref_hashSizesCount, $average) = &calculateStdDev($AverageSizeGenome, $lower_fence, $upper_fence, \@statistics);
	$FilterAverageGenome = $avg_filter;
	%SizesFrequency = %{$ref_hashSizesCount};
	$errorbarBegin = sprintf("%.2f",($AverageSizeGenome - $stddev)); $errorbarEnd = sprintf("%.2f",($AverageSizeGenome + $stddev));
	$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));

	if ( $type eq "introns" ) {
	    $content_introns = "";
	    $intronsCDSgenome = sprintf("%.2f",(($intronsCDS_sp{$sps} * 100) / $cds_sp));
	    $GeneDensityGenome = sprintf("%.2f",($number_sp{$sps}{$type} / $intronsCDS_sp{$sps}));

### weighted average length of introns in CDS bearing introns:
	    $MeanSizeGenome = sprintf("%.2f",($mean_sp{$sps}{$type} / $intronsCDS_sp{$sps}));
	    $IntronContentGenome = $GenomeFrequency;
            $content_introns = "($GenomeFrequencyMbs Mbs)";
	    $NfreqGenome_introns = $NfreqGenome;
            $ATfreqGenome_introns = $ATfreqGenome;
            $GCfreqGenome_introns = $GCfreqGenome;
	    $TranscriptStructure = sprintf("%.2f",($number_sp{$sps}{introns} / $number_sp{$sps}{CDSintrons}));	# $NumberBySpecies{"$genera $specie"}{CDSintrons}+= $numberExonsByCDS;
	    $IntronDensityExons = sprintf("%.2f",($mean_sp{$sps}{density} / $intronsCDS_sp{$sps}));		# $meanByCDSbySpecies{"$genera $specie"}{density}+= $densityIntronsByExons;
	}
	if ( $type eq "exons" ) {
	    $content_exons = "";
	    $GeneDensityGenome = sprintf("%.2f",($number_sp{$sps}{$type} / $cds_sp));

### weighted average length of exons in CDS:
	    $MeanSizeGenome = sprintf("%.2f",($mean_sp{$sps}{$type} / $cds_sp));
	    $ExonContentGenome = $GenomeFrequency;
	    $content_exons = "($GenomeFrequencyMbs Mbs)";
	    $NfreqGenome_exons = $NfreqGenome; $ATfreqGenome_exons = $ATfreqGenome; $GCfreqGenome_exons = $GCfreqGenome;
	}

	$logGenome = sprintf("%.4f", (log10($GenomeSizeMbs)));
	if ( $GenomeFrequencyMbs > 0 ) { $logSizeGenome = sprintf("%.4f", (log10($GenomeFrequencyMbs))); } else { $logSizeGenome = 0; }
	if ( $AverageSizeGenome > 0 ) { $logAverageGenome = sprintf("%.4f", (log10($AverageSizeGenome))); } else { $logAverageGenome = 0; }
	if ( $MeanSizeGenome > 0 ) { $logMeanGenome = sprintf("%.4f", (log10($MeanSizeGenome))); } else { $logMeanGenome = 0; }

	if ( $type eq "exons" || $type eq "introns" ) {
	    $WeightedAverageGenome = sprintf("%.2f",(($AverageSizeGenome + $MeanSizeGenome + $FilterAverageGenome) / 3 ));
	    @means = sort {$a <=> $b} @{ $wavg_sp{$sps}{$type} };
	    ($variance_wavg, $stddev_wavg, $varwavg_filter, $sdwavg_filter, $wavg_filter, $null, $waverage) = &calculateStdDev($MeanSizeGenome, $lower_fence, $upper_fence, \@means);
	    $errorbarBegin_wavg = sprintf("%.2f",($MeanSizeGenome - $stddev_wavg)); $errorbarEnd_wavg = sprintf("%.2f",($MeanSizeGenome + $stddev_wavg));
	    $errorbarBeginFilter_wavg = sprintf("%.2f",($wavg_filter - $sdwavg_filter)); $errorbarEndFilter_wavg = sprintf("%.2f",($wavg_filter + $sdwavg_filter));

	    @density = sort {$a <=> $b} @{ $density_sp{$sps}{$type} };
	    #print "$type = *@density*\n";
	    ($DQ1, $DQ3, $DQ_top, $DQ_down, $DIQR, $Dlower_fence, $Dupper_fence, $Dmedian) = &getQuartiles(\@density);
	    ($variance_wden, $stddev_wden, $varwden_filter, $sdwden_filter, $wavgden_filter, $ref_hashDensityCount, $daverage) = &calculateStdDev($GeneDensityGenome, $Dlower_fence, $Dupper_fence,\@density);
	    %DensityFrequency = %{$ref_hashDensityCount};
	    $errorbarBegin_wden = sprintf("%.2f",($GeneDensityGenome - $stddev_wden)); $errorbarEnd_wden = sprintf("%.2f",($GeneDensityGenome + $stddev_wden));
	    $errorbarBeginFilter_wden = sprintf("%.2f",($wavgden_filter - $sdwden_filter)); $errorbarEndFilter_wden = sprintf("%.2f",($wavgden_filter + $sdwden_filter));
	    $density_filter = sprintf("%.2f",($wavgden_filter));
	    #print "$type = density_total: $GeneDensityGenome = density_wave: $daverage = density_filter: $wavgden_filter\n";

	    @densitypos_filter = (); @density_filter = (); @densitypos_all = (); @density_all = (); @frequency = ();
	    foreach $sizes ( sort {$a <=> $b} keys %{ $DensityFrequency{filter} } ) {
	    	push @densitypos_filter, "$sizes"; push @density_filter, "$DensityFrequency{filter}{$sizes}";
		push @frequency, "$DensityFrequency{filter}{$sizes}";
		#print "*$sizes = $DensityFrequency{filter}{$sizes}*\t";
	    }
	    foreach $sizes ( sort {$a <=> $b} keys %{ $DensityFrequency{all} } ) {
	        push @densitypos_all, "$sizes"; push @density_all, "$DensityFrequency{all}{$sizes}";
		#print "*$sizes = $DensityFrequency{all}{$sizes}*\t";
	    }

	    if ( @density_all && @density_filter ) {
	        @avg_density = (); @avg_densityfilter = (); @weighted_density = (); @ratio_density = (); @high_density = (); @lower_ratio = ();
		$x_density = 0; $x_densityfilter = 0;
                $chart_density = ""; $x_all = ""; $x_filter = ""; $avg_density = ""; $avg_densityfilter = ""; $weighted_density = ""; $ratio_density = "";
		@frequency = sort {$a <=> $b} @frequency;
	        $x_density = $GeneDensityGenome - 1;
		$x_densityfilter = $wavgden_filter - 1;
	        push @avg_density, $x_density; push @high_density, $frequency[-1]; push @lower_ratio, 10;
		push @avg_densityfilter, $x_densityfilter; push @weighted_density, $IntronDensityExons; push @ratio_density, $TranscriptStructure;
		#print "lower ratio: @lower_ratio ------> ratio density: @ratio_density ------> weighted density: @weighted_density\n";

	        $chart_density = Chart::Gnuplot->new( 
				output  => "$outfile_gff.density.distribution.$type.genome.ps", 
				title   => { text => "Frequency of $type densities across the genome", color => "#000080", offset => "0, 2", }, 
				xlabel  => { text => "$type density (number of $type per protein-coding gene)", color => "#000080", offset => "0, -0.5", }, 
				ylabel  => { text => "frequency",  color => "#000080", }, 
				yrange  => [0, '*'], 
				xtics   => { font => ", 7", rotate => "90", mirror => 'off', },
				ytics   => { mirror => 'off', }, 
				y2range => [0, 1], 
				y2tics 	=> 'on',
				y2tics  => { mirror => 'off', }, 
				y2label => { text => "intron:exon ratio", color => "#000080", }, 
				legend  => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, },
				);
	        $x_all = Chart::Gnuplot::DataSet->new( xdata => \@densitypos_all, ydata => \@density_all, title => "all densities", color => "#7F7F7F", fill => { density => 0.8, border => "on", }, linetype => 'solid', style => "histograms", );
	        $x_filter = Chart::Gnuplot::DataSet->new( xdata => \@densitypos_filter, ydata => \@density_filter, title => "densities: lower-upper quartiles", color => "#000000", fill => { density => 0.8, border => "on", }, linetype => 'solid', style => "histograms", );
		$avg_density = Chart::Gnuplot::DataSet->new( xdata => \@avg_density, ydata => \@high_density, style => "impulses", linetype => "dot-longdash", color => "#0000FF", width => 2, title => "avg density ($GeneDensityGenome)", );
	     	$avg_densityfilter = Chart::Gnuplot::DataSet->new( xdata => \@avg_densityfilter, ydata => \@high_density, style => "impulses", linetype => "dot-longdash", color => "#FF1493", width => 2, title => "filtered avg density ($wavgden_filter)", );
		
		$weighted_density = Chart::Gnuplot::DataSet->new( 
				xdata => \@lower_ratio, ydata => \@weighted_density, style => "hlines", axes => "x2y2", linetype => "dot", color => "#32CD32", width => 2, title => "weighted intron:exon ratio ($IntronDensityExons)", 
				);
		$ratio_density = Chart::Gnuplot::DataSet->new( 
				xdata => \@lower_ratio, ydata => \@ratio_density, style => "hlines", axes => "x2y2", linetype => "dot", color => "#2F4F4F", width => 2, title => "exon:intron ratio ($TranscriptStructure)", 
				);

		if ( $type eq "exons" ) { $chart_density->plot2d($x_all, $x_filter, $avg_density, $avg_densityfilter); }
		if ( $type eq "introns" ) { $chart_density->plot2d($x_all, $x_filter, $avg_density, $avg_densityfilter, $weighted_density, $ratio_density); }
	    }
	}


### -------------------- ###
###  SECTION FOR EXONS
### -------------------- ###

        if ( $type eq "exons" ) {
	    @frequency = ();
            print STATSEXONS_GENOME "$sps\t$GenomeSizeMbs\t$cds_sp\t$number_sp{$sps}{$type}\t$GeneDensityGenome\t$density[0]\t$density[-1]\t$variance_wden\t$stddev_wden\t$errorbarBegin_wden\t$errorbarEnd_wden\t$density_filter\t$varwden_filter\t$sdwden_filter\t$errorbarBeginFilter_wden\t$errorbarEndFilter_wden\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
            #print "STATSEXONS_GENOME $sps\t$GenomeSizeMbs\t$cds_sp\t$number_sp{$sps}{$type}\t$GeneDensityGenome\t$density[0]\t$density[-1]\t$variance_wden\t$stddev_wden\t$errorbarBegin_wden\t$errorbarEnd_wden\t$density_filter\t$varwden_filter\t$sdwden_filter\t$errorbarBeginFilter_wden\t$errorbarEndFilter_wden\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";

	    &addQuartilesToCoordinates($outfile_gff, $type, $Q1, $Q3, $Q_top, $Q_down, $lower_fence, $upper_fence, $median, \%all_outfile);

	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{filter} } ) {
	        print SIZESEXONS_FILTER "$sizes\t$SizesFrequency{filter}{$sizes}\n";
	        push(@pairs_quartiles, [$sizes, $SizesFrequency{filter}{$sizes}]);
		push @frequency, "$SizesFrequency{filter}{$sizes}";
	    }
	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{all} } ) {
	        print SIZESEXONS_ALL "$sizes\t$SizesFrequency{all}{$sizes}\n";
	        push(@pairs_all, [$sizes, $SizesFrequency{all}{$sizes}]);
	    }
	    @frequency = sort {$a <=> $b} @frequency;
	}


### --------------------- ###
###  SECTION FOR INTRONS
### --------------------- ###

        if ( $type eq "introns" ) {
	    @frequency = (); $exon_shortintron = 0; $def_exons = 0; $def_total = 0;
            $deffilter_total = 0; $defexons_filter = 0; $exon_filterintron = 0;
	    &addQuartilesToCoordinates($outfile_gff, $type, $Q1, $Q3, $Q_top, $Q_down, $lower_fence, $upper_fence, $median, \%all_outfile);

	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{filter} } ) {
	        print SIZESINTRONS_FILTER "$sizes\t$SizesFrequency{filter}{$sizes}\n";
	        push(@pairs_quartiles, [$sizes, $SizesFrequency{filter}{$sizes}]);
		push @frequency, "$SizesFrequency{filter}{$sizes}";
	    }
	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{all} } ) {
	        print SIZESINTRONS_ALL "$sizes\t$SizesFrequency{all}{$sizes}\n";
	        push(@pairs_all, [$sizes, $SizesFrequency{all}{$sizes}]);
	    }
	    @frequency = sort {$a <=> $b} @frequency;

	    if ( @definition ) {
	        foreach $par ( 0 .. $#definition ) { 
		    push(@definition_all, ["$definition[$par][1]", "$definition[$par][2]"]);
		    $par1 = 0; $par2 = 0;
		    $par1 = $SizesFrequency{filter}{"$definition[$par][1]"} || 0;
		    $par2 = $SizesFrequency{filter}{"$definition[$par][2]"} || 0;
		    #print "$par1 ---> $par2\n";

		    if ( $par1 > 0  && $par2 > 0 ) {
		        #print "$definition[$par][0]\t*$definition[$par][1]*\t*$definition[$par][2]*\n";
		        print DEFINITION_FILTER "$definition[$par][0]\t$definition[$par][1]\t$definition[$par][2]\n";
		        push(@definition_filter, ["$definition[$par][1]", "$definition[$par][2]"]);
		    }
		    if ( $par1 > 0  || $par2 > 0  ) {
		        $deffilter_total ++;

			if ( $par1 <= 250 || $par2 <= 250 ) {
			    $defexons_filter ++;
			}
		    }
	        }
		#print "total filter: $deffilter_total || filtered: $defexons_filter\n";

		if ( $deffilter_total > 0 ) { $exon_filterintron = sprintf("%.2f",(($defexons_filter * 100) / $deffilter_total)); }
		($def_exons, $def_total) = split ":", $definition;
		if ( $def_total > 0 ) { $exon_shortintron = sprintf("%.2f",(($def_exons * 100) / $def_total)); }
		#print "$exon_shortintron -> $def_exons / $def_total -> $number_sp{$sps}{exons}\n";

	        $dataSet = ""; $chart3d_definition = "";
	        $chart3d_definition = Chart::Gnuplot->new( 
				output => "$outfile_gff.exons-introns.definition.sizes.ps", 
				title  => { text => "Distribution of lengths for all flanking introns of each exon in the genome", color => "#000080", offset => "0, 2", }, 
				xlabel => { text => "upstream intron size (nts)", color => "#000080", offset => "0, -1", }, 
				ylabel => { text => "downstream intron size (nts)", color => "#000080", }, 
				xtics  => { mirror => 'off', },
				ytics  => { mirror => 'off', },
				);
	        $dataSet = Chart::Gnuplot::DataSet->new( points => \@definition_all, style => "dots", width => 3, );
	        $chart3d_definition->plot2d($dataSet);

	        if ( @definition_filter ) {
	            $dataSet = ""; $chart3d_definition = ""; 
	            $chart3d_definition = Chart::Gnuplot->new( 
				output  => "$outfile_gff.quartiles.exons-introns.definition.sizes.ps", 
				title   => { text => "Distribution of lengths for filtered flanking introns of each exon in the genome", color => "#000080", }, 
				xlabel  => { text => "upstream intron size (nts)", color => "#000080", offset => "0, -1", }, 
				ylabel  => { text => "downstream intron size (nts)", color => "#000080", }, 
				x2label => { text => "$exon_shortintron% of exons ($def_total) are flanked by at least one short intron (<250 nts)", color => "#696969", offset => "0, 2", }, 
				xtics   => { mirror => 'off', },
				ytics   => { mirror => 'off', },
				);
	            $dataSet = Chart::Gnuplot::DataSet->new( points => \@definition_filter, style => "dots", width => 3, );
	            $chart3d_definition->plot2d($dataSet);
	        }
	    }
	    #$chart3d_definition->plot3d($dataSet);

            print STATSINTRONS_GENOME "$sps\t$GenomeSizeMbs\t$cds_sp\t$intronsCDSgenome\t$number_sp{$sps}{$type}\t$exon_shortintron\t$exon_filterintron\t$GeneDensityGenome\t$IntronDensityExons\t$TranscriptStructure\t$density[0]\t$density[-1]\t$variance_wden\t$stddev_wden\t$errorbarBegin_wden\t$errorbarEnd_wden\t$density_filter\t$varwden_filter\t$sdwden_filter\t$errorbarBeginFilter_wden\t$errorbarEndFilter_wden\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
            #print "STATSINTRONS_GENOME $sps\t$GenomeSizeMbs\t$cds_sp\t$intronsCDSgenome\t$number_sp{$sps}{$type}\t$exon_shortintron\t$exon_filterintron\t$GeneDensityGenome\t$IntronDensityExons\t$TranscriptStructure\t$density[0]\t$density[-1]\t$variance_wden\t$stddev_wden\t$errorbarBegin_wden\t$errorbarEnd_wden\t$density_filter\t$varwden_filter\t$sdwden_filter\t$errorbarBeginFilter_wden\t$errorbarEndFilter_wden\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";

	    for $exonNumber ( sort {$a <=> $b} keys %{ $order_exons{$sps} } ) {
	        @exonsizes = (); @filter = ();
		$sum = 0; $num = 0; $exonSizeByPosAve_filter = 0;
		$exonSizeByPositionAve = 0; $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $avg_filter = 0; $null = ""; 
		$Q_top = 0; $Q_down = 0; $Q1 = 0; $Q3 = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0; $median = 0; $average = 0;
		$errorbarBegin = 0; $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;

	        $exonSizeByPositionAve = sprintf("%.2f",( $order_exons{$sps}{$exonNumber}{summatory} / $order_exons{$sps}{$exonNumber}{number} ));
	        @exonsizes = sort {$a <=> $b} @{ $orderexons_stats{$sps}{$exonNumber} };	        
	        ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@exonsizes);
	        ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($exonSizeByPositionAve, $lower_fence, $upper_fence, \@exonsizes);
		$errorbarBegin = sprintf("%.2f",($exonSizeByPositionAve - $stddev)); $errorbarEnd = sprintf("%.2f",($exonSizeByPositionAve + $stddev));
		$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));
		#print "ORDER EXONS: *@exonsizes*\ntotal_ave: $exonSizeByPositionAve = check_total_ave: $average = filter_ave: $avg_filter\n";

	        push @exonpositions, "$exonNumber"; push @exonaverages, "$exonSizeByPositionAve"; push @exonerrors, "$stddev";
		print ORDER_EXONS "exon_$exonNumber\t$order_exons{$sps}{$exonNumber}{number}\t$exonSizeByPositionAve\t$avg_filter\t$exonsizes[0]\t$exonsizes[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\n";
		#print "ORDER_EXONS: exon_$exonNumber\t$order_exons{$sps}{$exonNumber}{number}\t$exonSizeByPositionAve\t$avg_filter\t$exonsizes[0]\t$exonsizes[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\n";

		foreach $filter ( @exonsizes ) {
		    if ( exists $SizesFrequency{filter}{$filter} ) { $sum += $filter; $num ++; push @filter, "$filter"; }
		}
		unless ( $num == 0 ) {
		    $errorbarBegin = 0; $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;
		    $exonSizeByPosAve_filter = sprintf("%.2f",( $sum / $num ));

		    $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $avg_filter = 0; $null = ""; 
		    $Q_top = 0; $Q_down = 0; $Q1 = 0; $Q3 = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0; $median = 0; $average = 0;
		    ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@filter);
	            ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($exonSizeByPosAve_filter, $lower_fence, $upper_fence, \@filter);
		    $errorbarBegin = sprintf("%.2f",($exonSizeByPosAve_filter - $stddev)); $errorbarEnd = sprintf("%.2f",($exonSizeByPosAve_filter + $stddev));
		    #print "ORDER FILTER EXONS: filter_total_ave: $exonSizeByPosAve_filter = filter_check_total_ave: $average = filter_filter_ave: $avg_filter\n";

		    push @exonpos_filter, "$exonNumber"; push @exonave_filter, "$exonSizeByPosAve_filter"; push @exonerrors_filter, "$stddev"; push @exon_counts, "$num";
		    print ORDEREXONS_FILTER "exon_$exonNumber\t$num\t$exonSizeByPosAve_filter\t$filter[0]\t$filter[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\n";
		    #print "ORDEREXONS_FILTER: exon_$exonNumber\t$num\t$exonSizeByPosAve_filter\t$filter[0]\t$filter[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\n";
		}
	    }
	    #print "exon positions: @exonpositions\nexon averages: @exonaverages\nexon errors: @exonerrors\n";
	    #print "exon positions filter: @exonpos_filter\nexon averages filter: @exonave_filter\n";

	    for $intronNumber ( sort {$a <=> $b} keys %{ $order_introns{$sps} } ) {
	    	@intronsizes = (); @filter = ();
		$null = ""; $sum = 0; $num = 0; $intronSizeByPosAve_filter = 0; $errorbarBegin = 0;
                $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;
		$intronSizeByPositionAve = 0; $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $avg_filter = 0;
		$Q_top = 0; $Q_down = 0; $Q1 = 0; $Q3 = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0; $median = 0; $average = 0;

		$intronSizeByPositionAve = sprintf("%.2f",( $order_introns{$sps}{$intronNumber}{summatory} / $order_introns{$sps}{$intronNumber}{number} ));
	        @intronsizes = sort {$a <=> $b} @{ $orderintrons_stats{$sps}{$intronNumber} };
	        ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@intronsizes);
	        ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($intronSizeByPositionAve, $lower_fence, $upper_fence, \@intronsizes);
		$errorbarBegin = sprintf("%.2f",($intronSizeByPositionAve - $stddev)); $errorbarEnd = sprintf("%.2f",($intronSizeByPositionAve + $stddev));
		$errorbarBeginFilter = sprintf("%.2f",($avg_filter - $sd_filter)); $errorbarEndFilter = sprintf("%.2f",($avg_filter + $sd_filter));
		#print "ORDER INTRONS: *@intronsizes*\ntotal_ave: $intronSizeByPositionAve = check_total_ave: $average = filter_ave: $avg_filter\n";

		push @intronpositions, "$intronNumber"; push @intronaverages, "$intronSizeByPositionAve"; push @intronerrors, "$stddev";
	        print ORDER_INTRONS "intron_$intronNumber\t$order_introns{$sps}{$intronNumber}{number}\t$intronSizeByPositionAve\t$avg_filter\t$intronsizes[0]\t$intronsizes[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\n";
	        #print "ORDER_INTRONS: intron_$intronNumber\t$order_introns{$sps}{$intronNumber}{number}\t$intronSizeByPositionAve\t$avg_filter\t$intronsizes[0]\t$intronsizes[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\n";

		foreach $filter ( @intronsizes ) {
		    if ( exists $SizesFrequency{filter}{$filter} ) { $sum += $filter; $num ++; push @filter, "$filter"; }
		}
		unless ( $num == 0 ) {
		    $errorbarBegin = 0; $errorbarEnd = 0; $errorbarBeginFilter = 0; $errorbarEndFilter = 0;
		    $intronSizeByPosAve_filter = sprintf("%.2f",( $sum / $num ));

		    $variance = 0; $stddev = 0; $var_filter = 0; $sd_filter = 0; $avg_filter = 0; $null = ""; 
		    $Q_top = 0; $Q_down = 0; $Q1 = 0; $Q3 = 0; $IQR = 0; $lower_fence = 0; $upper_fence = 0; $median = 0; $average = 0;
		    ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@filter);
	            ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null, $average) = &calculateStdDev($intronSizeByPosAve_filter, $lower_fence, $upper_fence, \@filter);
		    $errorbarBegin = sprintf("%.2f",($intronSizeByPosAve_filter - $stddev)); $errorbarEnd = sprintf("%.2f",($intronSizeByPosAve_filter + $stddev));
		    #print "ORDER FILTER INTRONS: filter_total_ave: $intronSizeByPosAve_filter = filter_check_total_ave: $average = filter_filter_ave: $avg_filter\n";

		    push @intronpos_filter, "$intronNumber"; push @intronave_filter, "$intronSizeByPosAve_filter"; push @intronerrors_filter, "$stddev"; push @intron_counts, "$num";
		    print ORDERINTRONS_FILTER "intron_$intronNumber\t$num\t$intronSizeByPosAve_filter\t$filter[0]\t$filter[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\n";
		    #print "ORDERINTRONS_FILTER: intron_$intronNumber\t$num\t$intronSizeByPosAve_filter\t$filter[0]\t$filter[-1]\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\n";
		}
	    }
	    #print "intron positions: @intronpositions\nintron averages: @intronaverages\nintron erros: @intronerrors\n";
	    #print "intron positions filter: @intronpos_filter\nintron averages filter: @intronave_filter\nintron count filter: @intron_counts\n";

	    ### PLOTS FOR EXON AND INTRON SIZES BY POSITION IN THE GENE (ALL POINTS AND ONLY QUARTILE POINTS)

	    $chart_order = ""; 	$x_exons = ""; $x_introns = ""; $scalar = "";
	    $scalar = scalar @intronpositions;
	    $chart_order = Chart::Gnuplot->new( 
				output    => "$outfile_gff.exons-introns.order.sizes.ps", 
				imagesize => "1.5, 1.0", 
				title     => { text => "Average lenght of all exons and introns by their position in a gene", color => "#000080", offset => "0, 2", }, 
				xlabel    => { text => "position of exon/intron within the gene", color => "#000080", offset => "0, -0.5", }, 
				ylabel    => { text => "average size (nts)", color => "#000080", }, 
				yrange    => [0, '*'], 
				xtics     => { font => ", 7", rotate => "90", mirror => 'off', },
				legend   => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, },
				);
	    $x_exons = Chart::Gnuplot::DataSet->new( 
				xdata => \@exonpositions, ydata => \@exonaverages, title => "exons", color => "#00008B", fill => { density => 0.7, border => "on", }, linetype => 'solid', style => "histograms", 
				);
	    $x_introns = Chart::Gnuplot::DataSet->new( 
				xdata => \@intronpositions, ydata => \@intronaverages, title => "introns", color => "red", fill => { density => 0.7, border => "on", }, linetype => 'solid', style => "histograms", 
				);
	    if ( @exonpositions && @intronpositions ) { $chart_order->plot2d($x_exons, $x_introns); }
	    elsif ( @exonpositions && $scalar == 0 ) { $chart_order->plot2d($x_exons); }

	    $chart_order = ""; 	$x_exons = ""; $x_introns = ""; $scalar = ""; $line_exons = ""; $line_introns = "";
            @exonpos_line = (); @intronpos_line = ();
	    $scalar = scalar @intronpos_filter;
	    @exonpos_line = @exonpos_filter; @intronpos_line = @intronpos_filter;
	    unshift @exonpos_line, "0"; pop @exonpos_line;
	    unshift @intronpos_line, "0"; pop @intronpos_line;
	    $chart_order = Chart::Gnuplot->new( 
				output    => "$outfile_gff.quartiles.exons-introns.order.sizes.ps", 
				imagesize => "1.5, 1.0", 
				title     => { text => "Average lenght of filtered exons and introns by their position in a gene", color => "#000080", offset => "0, 2", }, 
				xlabel    => { text => "position of exon/intron within the gene", color => "#000080", offset => "0, -0.5", }, 
				ylabel    => { text => "average size (nts)", color => "#000080", }, 
				yrange    => [0, '*'], 
				xtics     => { font => ", 7", rotate => "90", mirror => 'off', }, 
				ytics     => { mirror => 'off', },
				y2tics 	  => 'on',
				y2label   => { text => "exon/intron count", color => "#000080", }, 
				legend   => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, },
				);
	    $x_exons = Chart::Gnuplot::DataSet->new( 
				xdata => \@exonpos_filter, ydata => \@exonave_filter, title => "exons", color => "#00008B", fill => { density => 0.7, border => "on", }, linetype => 'solid', style => "histograms", 
				);
	    $x_introns = Chart::Gnuplot::DataSet->new( 
				xdata => \@intronpos_filter, ydata => \@intronave_filter, title => "introns", color => "red", fill => { density => 0.7, border => "on", }, linetype => 'solid', style => "histograms", 
				);
	    $line_exons = Chart::Gnuplot::DataSet->new( xdata => \@exonpos_line, ydata => \@exon_counts, title => "no. exons", style => "lines", axes => 'x1y2', color => "#00CED1", width => 8, );
	    $line_introns = Chart::Gnuplot::DataSet->new( xdata => \@intronpos_line, ydata => \@intron_counts, title => "no. introns", style => "lines", axes => 'x1y2', color => "#7FFF00", width => 8, );

	#-----------------
	    $chart_errors = ""; $exon_errors = ""; $intron_errors = "";
	    $chart_errors = Chart::Gnuplot->new( 
				output    => "$outfile_gff.errorbars.quartiles.exons-introns.order.sizes.ps", 
				imagesize => "1.5, 1.0", 
				title     => { text => "Average lenght of filtered exons and introns by their position in a gene", color => "#000080", offset => "0, 2", }, 
				xlabel    => { text => "position of exon/intron within the gene", color => "#000080", offset => "0, -0.5", }, 
				ylabel    => { text => "average size (nts)", color => "#000080", }, 
				yrange    => [0, '*'], 
				xtics     => { rotate => "90", mirror => 'off', },
				);
	    $exon_errors = Chart::Gnuplot::DataSet->new(
	    	    		xdata     => \@exonpos_filter, ydata => [[@exonave_filter], [@exonerrors_filter]], title => "exons", color => "#00008B", width => 3, style => "yerrorbars", 
	    			);
	    $intron_errors = Chart::Gnuplot::DataSet->new(
	    			xdata     => \@intronpos_filter, ydata => [[@intronave_filter], [@intronerrors_filter]], title => "introns", color => "red", width => 3, style => "yerrorbars", 
				);

	    if ( @exonpos_filter && @intronpos_filter ) { $chart_order->plot2d($x_exons, $x_introns, $line_exons, $line_introns); $chart_errors->plot2d($exon_errors, $intron_errors); }
	    elsif ( @exonpos_filter && $scalar == 0 ) { $chart_order->plot2d($x_exons, $line_exons); $chart_errors->plot2d($exon_errors); }
	}

### ------------------------- ###
###  SECTION FOR INTERGENICS
### ------------------------- ###

        if ( $type eq "intergenics" ) {
	    @frequency = (); $content_intergenics = "";
	    $MeanSizeGenome = $avg_filter;
	    $IntergenicContentGenome = $GenomeFrequency; $content_intergenics = "($GenomeFrequencyMbs Mbs)";
	    $NfreqGenome_intergenics = $NfreqGenome; $ATfreqGenome_intergenics = $ATfreqGenome; $GCfreqGenome_intergenics = $GCfreqGenome;
	    $WeightedAverageGenome = sprintf("%.2f",(($AverageSizeGenome + $FilterAverageGenome) / 2 ));

	    print STATSINTERGENICS_GENOME "$sps\t$GenomeSizeMbs\t$cds_sp\t$GenomeFrequencyMbsCDS\t$GenomeFrequencyCDS\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\n";
	    #print "STATSINTERGENICS_GENOME $sps\t$GenomeSizeMbs\t$cds_sp\t$GenomeFrequencyMbsCDS\t$GenomeFrequencyCDS\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$length_sp{$sps}{$type}\t$GenomeFrequencyMbs\t$GenomeFrequencyTotal\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$FilterAverageGenome\t$WeightedAverageGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\n";

	    &addQuartilesToCoordinates($outfile_gff, $type, $Q1, $Q3, $Q_top, $Q_down, $lower_fence, $upper_fence, $median, \%all_outfile);

	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{filter} } ) {
	    	print SIZESINTERGENICS_FILTER "$sizes\t$SizesFrequency{filter}{$sizes}\n";
		push(@pairs_quartiles, [$sizes, $SizesFrequency{filter}{$sizes}]);
		push @frequency, "$SizesFrequency{filter}{$sizes}";
	    }
	    foreach $sizes ( sort {$a <=> $b} keys %{ $SizesFrequency{all} } ) {
	    	print SIZESINTERGENICS_ALL "$sizes\t$SizesFrequency{all}{$sizes}\n";
		push(@pairs_all, [$sizes, $SizesFrequency{all}{$sizes}]);
	    }
	    @frequency = sort {$a <=> $b} @frequency;

            ### -------------------- ###
            ### PLOTS FOR EXON - INTRON - INTERGENIC LENGTH RANGES
            ### -------------------- ###

	    $chart_ranges = ""; $r10 = ""; $r50 = ""; $r100 = ""; $r250 = ""; $r500 = ""; $r1000 = ""; $r5000 = ""; $r10000 = ""; $r50000 = ""; $r100000 = ""; $r500000 = "";
	    $r1000000 = ""; $r5000000 = ""; $r10000000 = ""; $r50000000 = ""; $r100000000 = ""; $r500000000 = ""; $r1000000000 = ""; 
	    @x_ranges = ( "exons_{$sum_exons}", "introns_{$sum_introns}", "intergenics_{$sum_intergenics}" );

	    @y_10 = (); @y_50 = (); @y_100 = (); @y_250 = (); @y_500 = (); @y_1000 = (); @y_5000 = (); @y_10000 = (); @y_50000 = (); @y_100000 = (); @y_500000 = ();
	    @y_1000000 = (); @y_5000000 = (); @y_10000000 = (); @y_50000000 = ();  @y_100000000 = (); @y_500000000 = (); @y_1000000000 = ();


    	    @y_10 = ( "$y_exons[0]", "$y_introns[0]", "$y_intergenics[0]" );
	    @y_50 = ( "$y_exons[1]", "$y_introns[1]", "$y_intergenics[1]" );
	    @y_100 = ( "$y_exons[2]", "$y_introns[2]", "$y_intergenics[2]" );
	    @y_250 = ( "$y_exons[3]", "$y_introns[3]", "$y_intergenics[3]" );
	    @y_500 = ( "$y_exons[4]", "$y_introns[4]", "$y_intergenics[4]" );
	    @y_1000 = ( "$y_exons[5]", "$y_introns[5]", "$y_intergenics[5]" );
	    @y_5000 = ( "$y_exons[6]", "$y_introns[6]", "$y_intergenics[6]" );
	    @y_10000 = ( "$y_exons[7]", "$y_introns[7]", "$y_intergenics[7]" );
	    @y_50000 = ( "$y_exons[8]", "$y_introns[8]", "$y_intergenics[8]" );
	    @y_100000 = ( "$y_exons[9]", "$y_introns[9]", "$y_intergenics[9]" );
	    @y_500000 = ( "$y_exons[10]", "$y_introns[10]", "$y_intergenics[10]" );
	    @y_1000000 = ( "$y_exons[11]", "$y_introns[11]", "$y_intergenics[11]" );
	    @y_5000000 = ( "$y_exons[12]", "$y_introns[12]", "$y_intergenics[12]" );
	    @y_10000000 = ( "$y_exons[13]", "$y_introns[13]", "$y_intergenics[13]" );
	    @y_50000000 = ( "$y_exons[14]", "$y_introns[14]", "$y_intergenics[14]" );
	    @y_100000000 = ( "$y_exons[15]", "$y_introns[15]", "$y_intergenics[15]" );
	    @y_500000000 = ( "$y_exons[16]", "$y_introns[16]", "$y_intergenics[16]" );
	    @y_1000000000 = ( "$y_exons[17]", "$y_introns[17]", "$y_intergenics[17]" );
	    #print "@x_ranges\n@y_10\n@y_50\n@y_100\n@y_250\n@y_500\n@y_1000\n@y_5000\n@y_10000\n@y_50000\n@y_100000\n@y_500000\n@y_1000000\n@y_5000000\n@y_10000000\n@y_50000000\n@y_100000000\n@y_500000000\n@y_1000000000\n";

    	    $chart_ranges = Chart::Gnuplot->new( 
	    			output   => "$output.ranges.sizes.distribution.genome.ps", 
				title    => { text => "Distribution of nucleotide ranges from the genomic features across the genome", color => "#000080", offset => "0, 2", }, 
				yrange   => [0, 100], 
				ylabel 	 => { text => "percentage from the total genomic features", color => "#000080", }, 
				boxwidth => "0.6 relative", 
				offset 	 => "-0.6,-0.6,0,0", 
				legend   => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, },
				"style histogram" => "rowstacked", 
				);
      	    #$r10 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_10, title => "0-10", color => "#000000", fill => {density => 0.9}, style => "histograms", );
    	    $r50 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_50, title => "15-50",  color => "#8B4513", fill => {density => 0.9}, style => "histograms", );
    	    $r100 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_100, title => "51-100", color => "#FFD700", fill => {density => 0.9}, style => "histograms", );
    	    $r250 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_250, title => "101-250", color => "#008000", fill => {density => 0.9}, style => "histograms", );
    	    $r500 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_500, title => "251-500", color => "#8B0000", fill => {density => 0.9}, style => "histograms", );
    	    $r1000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_1000, title => "501-1000", color => "#00CED1", fill => {density => 0.9}, style => "histograms", );
    	    $r5000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_5000, title => "1001-5000", color => "#00008B", fill => {density => 0.9}, style => "histograms", );
    	    $r10000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_10000, title => "10^{3}", color => "#C71585", fill => {density => 0.9}, style => "histograms", );
    	    $r50000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_50000, title => "50^{3}", color => "#808000", fill => {density => 0.9}, style => "histograms", );
    	    $r100000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_100000, title => "10^{4}", color => "#FFC0CB", fill => {density => 0.9}, style => "histograms", );
    	    $r500000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_500000, title => "50^{4}", color => "#8B0000", fill => {density => 0.9}, style => "histograms", );
    	    $r1000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_1000000, title => "10^{5}", color => "#ADFF2F", fill => {density => 0.9}, style => "histograms", );
    	    $r5000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_5000000, title => "50^{5}", color => "#00FFFF", fill => {density => 0.9}, style => "histograms", );
    	    $r10000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_10000000, title => "10^{6}", color => "#FF4500", fill => {density => 0.9}, style => "histograms", );
    	    $r50000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_50000000, title => "50^{6}", color => "#2F4F4F", fill => {density => 0.9}, style => "histograms", );
    	    $r100000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_100000000, title => "10^{7}", color => "#FF00FF", fill => {density => 0.9}, style => "histograms", );
    	    $r500000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_500000000, title => "50^{7}", color => "#8A2BE2", fill => {density => 0.9}, style => "histograms", );
    	    $r1000000000 = Chart::Gnuplot::DataSet->new( xdata => \@x_ranges, ydata => \@y_1000000000, title => "10^{8}", color => "#008B8B", fill => {density => 0.9}, style => "histograms", );
	    $chart_ranges->plot2d($r50, $r100, $r250, $r500, $r1000, $r5000, $r10000, $r50000, $r100000, $r500000, $r1000000, $r5000000, $r10000000, $r50000000, $r100000000, $r500000000, $r1000000000);
	    #timestamp => { fmt => '%d/%m/%y %H:%M', },


            ### -------------------- ###
            ### PLOTS FOR EXON - INTRON - INTERGENIC CONTENT AND NUCLEOTIDES
            ### -------------------- ###

	    my @x_features = (); my @content = (); my @cg_genome = (); my @at_genome = (); my @n_genome = ();
	    my $chart_content = ""; my $content = ""; my $cg_genome = ""; my $at_genome = ""; my $n_genome = "";

	    @x_features = ("exons_{$content_exons}", "introns_{$content_introns}", "intergenics_{$content_intergenics}");
	    #print "$ExonContentGenome ---> $IntronContentGenome ---> $IntergenicContentGenome\n";

	    @content = ($ExonContentGenome, $IntronContentGenome, $IntergenicContentGenome);
	    @cg_genome = ($GCfreqGenome_exons, $GCfreqGenome_introns, $GCfreqGenome_intergenics);
	    @at_genome = ($ATfreqGenome_exons, $ATfreqGenome_introns, $ATfreqGenome_intergenics);
	    @n_genome = ($NfreqGenome_exons, $NfreqGenome_introns, $NfreqGenome_intergenics);

    	    $chart_content = Chart::Gnuplot->new( 
	    			output  => "$output.features.content.distribution.genome.ps", 
				title   => { text => "Percentage of genome content for every feature in $sps", color => "#000080", offset => "0, 2", }, 
				yrange  => [0, 100], 
				ylabel  => { text => "percentage from the total genome size ($GenomeSizeMbs Mbs)", color => "#000080", }, 
				offset  => "-0.6,-0.6,0,0", 
				legend  => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, }, 
				y2tics 	=> 'on',
				y2label => { text => "percentage from the total feature size", color => "#000080", }, 
				xtics   => { mirror => 'off', },
				);
	    $content = Chart::Gnuplot::DataSet->new( xdata => \@x_features, ydata => \@content, title => "exon/intron/intergenic content", color => "#800080", fill => {density => 0.9}, style => "histograms", );
	    $cg_genome = Chart::Gnuplot::DataSet->new( xdata => \@x_features, ydata => \@cg_genome, title => "GC feature content", color => "#2F4F4F", fill => {density => 0.9}, style => "histograms", );
	    $at_genome = Chart::Gnuplot::DataSet->new( xdata => \@x_features, ydata => \@at_genome, title => "AT feature content", color => "#4682B4", fill => {density => 0.9}, style => "histograms", );
	    $n_genome = Chart::Gnuplot::DataSet->new( xdata => \@x_features, ydata => \@n_genome, title => "Ns feature content", color => "#7F7F7F", fill => {density => 0.9}, style => "histograms", );
	    $chart_content->plot2d($content, $cg_genome, $at_genome, $n_genome);
	}


        ### -------------------- ###
        ### PLOTS FOR EXON AND INTRON SIZES (ALL POINTS AND ONLY QUARTILE POINTS)
        ### -------------------- ###

	if ( @pairs_all ) {
	     my $chart_sizes = ""; my $dot_sizes = "";
	     $chart_sizes = Chart::Gnuplot->new( 
				output => "$output.all.sizes.$type.distribution.genome.ps", 
				title  => { text => "Frequency of all $type sizes across the genome", color => "#000080", offset => "0, 2", }, 
				xlabel => { text => "$type size (nts)", color => "#000080", offset => "0, -1", }, 
    				ylabel => { text => "frequency", color => "#000080", }, 
				xtics  => { mirror => 'off', },
				ytics  => { mirror => 'off', },
				);
	     $dot_sizes = Chart::Gnuplot::DataSet->new( points => \@pairs_all, style  => "dots", width => 3, );
	     $chart_sizes->plot2d($dot_sizes);
	}
	if ( @pairs_quartiles ) {
	     my @avg_size = (); my @mean_size = (); my @favg_size = (); my @wavg_size = (); my @high_avg = (); my @high_mean = (); my @high_favg = (); my @high_wavg = ();
	     my $avg_size = ""; my $mean_size = ""; my $favg_size = ""; my $wavg_size = ""; my $chart_sizes = ""; my $dot_sizes = ""; $chart_sizes = ""; $dot_sizes = "";

	     push @avg_size, $AverageSizeGenome; push @high_avg, $frequency[-1];
	     push @mean_size, $MeanSizeGenome; push @high_mean, $frequency[-1];
	     push @favg_size, $FilterAverageGenome; push @high_favg, $frequency[-1];
	     push @wavg_size, $WeightedAverageGenome; push @high_wavg, $frequency[-1];

	     $chart_sizes = Chart::Gnuplot->new( 
				output => "$output.quartiles.sizes.$type.distribution.genome.ps", 
				title  => { text => "Frequency of $type sizes within the lower and upper quartiles across the genome", color => "#000080", offset => "0, 2", }, 
				xlabel => { text => "$type size (nts)", color => "#000080", offset => "0, -1", }, 
    				ylabel => { text => "frequency", color => "#000080", }, 
				xtics  => { mirror => 'off', },
				ytics  => { mirror => 'off', },
				legend => { position => "outside center bottom", order => "horizontal reverse", align => "right", key => { font => ",3", }, },
				);
	     $dot_sizes = Chart::Gnuplot::DataSet->new( points => \@pairs_quartiles, style  => "dots", width => 3, );
	     $avg_size  = Chart::Gnuplot::DataSet->new( xdata => \@avg_size, ydata => \@high_avg, style => "impulses", linetype => "solid", color => "#2F4F4F", width => 2, title => "avg size ($AverageSizeGenome)", );
	     $mean_size = Chart::Gnuplot::DataSet->new( xdata => \@mean_size, ydata => \@high_mean, style => "impulses", linetype => "solid", color => "#800080", width => 2, title => "weighted avg size ($MeanSizeGenome)", );
	     $favg_size = Chart::Gnuplot::DataSet->new( xdata => \@favg_size, ydata => \@high_favg, style => "impulses", linetype => "solid", color => "#32CD32", width => 2, title => "quartile avg size ($FilterAverageGenome)", );
	     $wavg_size = Chart::Gnuplot::DataSet->new( xdata => \@wavg_size, ydata => \@high_wavg, style => "impulses", linetype => "solid", color => "#0000CD", width => 2, title => "normalised avg size ($WeightedAverageGenome)", );
	     if ( $type eq "intergenics" ) { $chart_sizes->plot2d($dot_sizes, $avg_size, $favg_size, $wavg_size); }
	     else { $chart_sizes->plot2d($dot_sizes, $avg_size, $mean_size, $favg_size, $wavg_size); }
	}

### -------------------- ###
### PRINTING GLOBAL DESCRIPTORS
### -------------------- ###

        if ( $sequence eq "assembly" && $type eq "exons" ) {
            print STATSEXONS_CHROMOSOME "$sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$GeneDensityGenome\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
            #print "STATSEXONS_CHROMOSOME: $sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$GeneDensityGenome\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$var_filter\t$sd_filter\t$variance_wavg\t$stddev_wavg\t$varwavg_filter\t$sdwavg_filter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
        }
        if ( $sequence eq "assembly" && $type eq "introns" ) {
            print STATSINTRONS_CHROMOSOME "$sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$intronsCDSgenome\t$GeneDensityGenome\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$variance_wavg\t$stddev_wavg\t$errorbarBegin_wavg\t$errorbarEnd_wavg\t$varwavg_filter\t$sdwavg_filter\t$errorbarBeginFilter_wavg\t$errorbarEndFilter_wavg\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
            #print "STATSINTRONS_CHROMOSOME: $sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$intronsCDSgenome\t$GeneDensityGenome\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$MeanSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$var_filter\t$sd_filter\t$variance_wavg\t$stddev_wavg\t$varwavg_filter\t$sdwavg_filter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\t$logMeanGenome\n";
        }
        if ( $sequence eq "assembly" && $type eq "intergenics" ) {
    	    print STATSINTERGENICS_CHROMOSOME "$sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$GenomeFrequencyMbsCDS\t$GenomeFrequencyCDS\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$errorbarBegin\t$errorbarEnd\t$var_filter\t$sd_filter\t$errorbarBeginFilter\t$errorbarEndFilter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\t$logAverageGenome\n";
    	    #print "STATSINTERGENICS_CHROMOSOME: $sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t$GenomeFrequencyMbsCDS\t$GenomeFrequencyCDS\t$number_sp{$sps}{$type}\t$LengthMbsGenome\t$GenomeFrequencyMbs\t$GenomeFrequency\t$statistics[0]\t$statistics[-1]\t$AverageSizeGenome\t$Q_down\t$Q1\t$median\t$Q3\t$Q_top\t$lower_fence\t$upper_fence\t$variance\t$stddev\t$var_filter\t$sd_filter\t$AfreqGenome\t$TfreqGenome\t$GfreqGenome\t$CfreqGenome\t$NfreqGenome\t$ATfreqGenome\t$GCfreqGenome\t$logGenome\t$logSizeGenome\n";
        }
    }

    if ( $type eq "introns" && !defined $number_sp{$sps}{introns} ) {
        print STATSINTRONS_GENOME "$sps\t$GenomeSizeMbs\t$cds_sp\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
        #print "STATSINTRONS_GENOME: $sps\t$GenomeSizeMbs\t$cds_sp\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";

        print STATSINTRONS_CHROMOSOME "$sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
        #print "STATSINTRONS_CHROMOSOME: $sps\tassembly\t$GenomeSizeMbs\t$cds_sp\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
    }
    print "**$type -> genome finished**\n";
}




######################################################################################################
#                                                                                                    #
#   CALCULATING SEVERAL STATISTIC INDIVIDUAL DESCRIPTORS                                             #
#                                                                                                    #
#   - freqATGC        ---> calculate the frequency of nucleotides in a sequence                      #
#   - getIntronType   ---> estimate the intron type of length modulo 3                               #
#   - getFreqATGCtype ---> calculate the fraction of nucleotide types from the total sequence size   #
#                                                                                                    # 
######################################################################################################


sub freqATGC {
    # INPUT: ($ref_hashByATGCgenome, $total_nts) = &freqATGC(\$$FastaByChromosomes{$seqname});
    my %hashByATGC = (); my $atgc= ""; my $nts = 0; my $total = 0;
    $atgc = shift;

    $hashByATGC{A} = 0; $hashByATGC{T} = 0; $hashByATGC{G} = 0; $hashByATGC{C} = 0; $hashByATGC{N} = 0;

    $total = length($$atgc);

    $hashByATGC{A} = ($$atgc =~ tr/A//);
    $hashByATGC{T} = ($$atgc =~ tr/T//);
    $hashByATGC{G} = ($$atgc =~ tr/G//);
    $hashByATGC{C} = ($$atgc =~ tr/C//);
    $hashByATGC{N} = ($$atgc =~ tr/N//);

    #$nts = $hashByATGC{A} + $hashByATGC{T} + $hashByATGC{G} + $hashByATGC{C} + $hashByATGC{N};
    #print "total = $total ---> nts = $nts\n";

    return (\%hashByATGC, $total);
    $atgc= ""; $nts = 0; $total = 0; %hashByATGC = ();
}


sub getIntronType {
    # INPUT: $multipleIntron = &getIntronType($intronSize);

    my ($size, $multiple, $intron_3n, $intron_type, $number) = (0, 0, 0, "", 0);
    $size =  shift;

    if ( $size >= 3 ) {
        $multiple = int($size / 3);
        $intron_3n =  $multiple * 3;
	$number = $size - $intron_3n;
        if ( $number == 0 ) { $intron_type = "3n"; return $intron_type; }
        if ( $number == 1 ) { $intron_type = "3n+1"; return $intron_type; }
        if ( $number == 2 ) { $intron_type = "3n+2"; return $intron_type; }
    }
    else { return $size;}
    ($size, $multiple, $intron_3n, $intron_type) = (0, 0, 0, "");
}


sub getFreqATGCtype {
    # INPUT: getFreqATGCtype("chromosomes", $species, $seqname, $type, \%atgc_chr);

    my %atgc = ();
    my $option = ""; my $name = ""; my $fastaname = ""; my $class = ""; my $plusORminus = "";
    my $TotalNTS = 0; my $Afreq = 0; my $Tfreq = 0; my $Gfreq = 0; my $Cfreq = 0; my $Nfreq = 0; my $ATfreq = 0; my $GCfreq = 0;

    $option = shift; $name = shift; $fastaname = shift; $class = shift;
    %atgc = %{$_[0]};

    if ( $option eq "chromosomes" ) {
        $TotalNTS = $atgc{$name}{$fastaname}{$class}{A} + $atgc{$name}{$fastaname}{$class}{T} + $atgc{$name}{$fastaname}{$class}{G} + $atgc{$name}{$fastaname}{$class}{C} + $atgc{$name}{$fastaname}{$class}{N};
    	$Afreq = sprintf("%.2f",(($atgc{$name}{$fastaname}{$class}{A} * 100) / $TotalNTS));
    	$Tfreq = sprintf("%.2f",(($atgc{$name}{$fastaname}{$class}{T} * 100) / $TotalNTS));
    	$Gfreq = sprintf("%.2f",(($atgc{$name}{$fastaname}{$class}{G} * 100) / $TotalNTS));
    	$Cfreq = sprintf("%.2f",(($atgc{$name}{$fastaname}{$class}{C} * 100) / $TotalNTS));
    	$Nfreq = sprintf("%.2f",(($atgc{$name}{$fastaname}{$class}{N} * 100) / $TotalNTS));
    	$ATfreq = sprintf("%.2f",((($atgc{$name}{$fastaname}{$class}{A} + $atgc{$name}{$fastaname}{$class}{T}) * 100) / $TotalNTS));
    	$GCfreq = sprintf("%.2f",((($atgc{$name}{$fastaname}{$class}{G} + $atgc{$name}{$fastaname}{$class}{C}) * 100) / $TotalNTS));

	return ($TotalNTS, $Afreq, $Tfreq, $Gfreq, $Cfreq, $Nfreq, $ATfreq, $GCfreq);
    }

    if ( $option eq "+" || $option eq "-" ) {
        $plusORminus = $option;
    	$TotalNTS = $atgc{$name}{$class}{$plusORminus}{A} + $atgc{$name}{$class}{$plusORminus}{T} + $atgc{$name}{$class}{$plusORminus}{G} + $atgc{$name}{$class}{$plusORminus}{C} + $atgc{$name}{$class}{$plusORminus}{N};
	#print "strands: *$plusORminus* -$TotalNTS-\n";
    	$Afreq = sprintf("%.2f",(($atgc{$name}{$class}{$plusORminus}{A} * 100) / $TotalNTS));
    	$Tfreq = sprintf("%.2f",(($atgc{$name}{$class}{$plusORminus}{T} * 100) / $TotalNTS));
    	$Gfreq = sprintf("%.2f",(($atgc{$name}{$class}{$plusORminus}{G} * 100) / $TotalNTS));
    	$Cfreq = sprintf("%.2f",(($atgc{$name}{$class}{$plusORminus}{C} * 100) / $TotalNTS));
    	$Nfreq = sprintf("%.2f",(($atgc{$name}{$class}{$plusORminus}{N} * 100) / $TotalNTS));
    	$ATfreq = sprintf("%.2f",((($atgc{$name}{$class}{$plusORminus}{A} + $atgc{$name}{$class}{$plusORminus}{T}) * 100) / $TotalNTS));
    	$GCfreq = sprintf("%.2f",((($atgc{$name}{$class}{$plusORminus}{G} + $atgc{$name}{$class}{$plusORminus}{C}) * 100) / $TotalNTS));

	return ($TotalNTS, $Afreq, $Tfreq, $Gfreq, $Cfreq, $Nfreq, $ATfreq, $GCfreq);
    }

    if ( $option eq "genome" ) {
        $TotalNTS = $atgc{$name}{$class}{A} + $atgc{$name}{$class}{T} + $atgc{$name}{$class}{G} + $atgc{$name}{$class}{C} + $atgc{$name}{$class}{N};
        $Afreq = sprintf("%.2f",(($atgc{$name}{$class}{A} * 100) / $TotalNTS));
        $Tfreq = sprintf("%.2f",(($atgc{$name}{$class}{T} * 100) / $TotalNTS));
	$Gfreq = sprintf("%.2f",(($atgc{$name}{$class}{G} * 100) / $TotalNTS));
        $Cfreq = sprintf("%.2f",(($atgc{$name}{$class}{C} * 100) / $TotalNTS));
        $Nfreq = sprintf("%.2f",(($atgc{$name}{$class}{N} * 100) / $TotalNTS));
        $ATfreq = sprintf("%.2f",((($atgc{$name}{$class}{A} + $atgc{$name}{$class}{T}) * 100) / $TotalNTS));
        $GCfreq = sprintf("%.2f",((($atgc{$name}{$class}{G} + $atgc{$name}{$class}{C}) * 100) / $TotalNTS));

	return ($TotalNTS, $Afreq, $Tfreq, $Gfreq, $Cfreq, $Nfreq, $ATfreq, $GCfreq);
    }
    %atgc = (); $option = ""; $name = ""; $fastaname = ""; $class = ""; $plusORminus = "";
    $TotalNTS = 0; $Afreq = 0; $Tfreq = 0; $Gfreq = 0; $Cfreq = 0; $Nfreq = 0; $ATfreq = 0; $GCfreq = 0;
}




###########################################################################################
#                                                                                         #
#   MAPPING REFERENCE FEATURE SETS ONTO GENOME SEQUENCE &                                 #
#   ESTIMATING GENOME COVERAGE FOR EACH FEATURE SET                                       #
#                                                                                         #
#   - getOverlaps ---> It calculates the net nucleotide lenght	                          #
#                      of the corresponding union intervals of:                           #
#                      1) all exons (EpcU)                                                #
#                      2) all introns flanked by protein-coding exons (IpcU)              #
#                      3) all intergenic regions (IRU)                                    #
#                                                                                         #
#   - getIntersections ---> It calculates the final length of introns:                    #
#                           the IpcU intervals that do not overlap exonic nucleotides     #
#                                                                                         #
#   - getNonCDSoverlaps ---> It calculates the net nucleotide lenght                      #
#                            of all union intervals of intergenic regions                 #
#                            (IRU) that do not overlap the exonic and intronic contents   #
#                                                                                         #
###########################################################################################


sub getOverlaps {
    # INPUT & OUTPUT: ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getOverlaps("$size_chr{$sps}{$seqname}", \@coordinates);

    my $positions= (); my @genomeOverlap = (); my @overlaps = (); my @lengths = ();
    my $pair = ""; my $leftpos = 0; my $rightpos = 0; my $size = 0;
    my $i = 0; my $z = 0;

    $size = shift;
    $positions = shift;

    for ( $z = 0; $z < ($size+1); $z++ ) {
	push @genomeOverlap, "0";
    }

    foreach $pair ( @$positions ) {
        chomp $pair;
	$leftpos = ""; $rightpos = "";
        ($leftpos, $rightpos) = split ":", $pair;
        for ( $i = $leftpos; $i <= $rightpos; $i++ ) { $genomeOverlap[$i] = $genomeOverlap[$i] + 1; }
    }
    $positions= ();
    my $SizecodingOverlap = 0;
    my $Sizenocoding = 0;
    my $total = 0;
    my $start = 0; my $end = 0;
    my $x = 0; my $y = 0;

    for ( $i = 1; $i <= $#genomeOverlap; $i++ ) {
        if ( $genomeOverlap[$i] > 0 && $SizecodingOverlap == 0 ) {
            $start = $i;
            $SizecodingOverlap++;
            $Sizenocoding = 0;
        }
	elsif ( $genomeOverlap[$i] == 0 && $Sizenocoding == 0 ) {
            if ( $SizecodingOverlap > 0 && $i < $#genomeOverlap ) {
		$end = $i - 1;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
            }
            $SizecodingOverlap = 0;
            $Sizenocoding++;
        }
	elsif ( $SizecodingOverlap > 0 && $i == $#genomeOverlap ) {
	        $end = $i;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
        }
    }
    return ($total, \@overlaps, \@lengths);
    @genomeOverlap = (); @overlaps = (); @lengths = (); $total = 0;
}


sub getIntersections {
    # INPUT: ($totalSizeIntron, $ref_arrayOverlapsIntrons, $ref_arrayLengthIntrons) = &getIntersections("$size_chr{$sps}{$seqname}", \@OverlapsIntrons, \@OverlapsExons);

    my $introns = (); my $exons = (); my @genomeOverlap = (); my @overlaps = (); my @lengths = ();
    my $pair = ""; my $leftpos = 0; my $rightpos = 0; my $size = 0;
    my $i = 0; my $z = 0;

    $size = shift;
    $introns = shift; $exons = shift;

    for ( $z = 0; $z < ($size+1); $z++ ) {
	push @genomeOverlap, "0";
    }

    #print "INTRONS: \t";
    foreach $pair ( @$introns ) {
        chomp $pair;
	#print "$pair\t";
	$leftpos = ""; $rightpos = "";
        ($leftpos, $rightpos) = split ":", $pair;
        for ( $i = $leftpos; $i <= $rightpos; $i++ ) { $genomeOverlap[$i] = $genomeOverlap[$i] + 1; }
    }

    #print "\nEXONS\t";
    foreach $pair ( @$exons ) {
        chomp $pair;
	#print "$pair\t";
	$leftpos = ""; $rightpos = "";
        ($leftpos, $rightpos) = split ":", $pair;
        for ( $i = $leftpos; $i <= $rightpos; $i++ ) { $genomeOverlap[$i] = "X"; }
    }
    #print "\n";
    #print "@genomeOverlap\n";
    $introns = (); $exons = (); 

    my $SizecodingOverlap = 0;
    my $Sizenocoding = 0;
    my $total = 0;
    my $start = 0; my $end = 0;
    my $x = 0; my $y = 0;

    for ( $i = 1; $i <= $#genomeOverlap; $i++ ) {
        if ( ($genomeOverlap[$i] ne "X") && ($genomeOverlap[$i] > 0) && ($SizecodingOverlap == 0) ) {
            $start = $i;
            $SizecodingOverlap++;
            $Sizenocoding = 0;
        }
	elsif ( ($genomeOverlap[$i] eq "X") || ($genomeOverlap[$i] == 0) && ($Sizenocoding == 0) ) {
            if ( $SizecodingOverlap > 0 && $i < $#genomeOverlap ) {
	        $end = $i - 1;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
            }
            $SizecodingOverlap = 0;
            $Sizenocoding++;
        }
	elsif ( $SizecodingOverlap > 0 && $i == $#genomeOverlap ) {
	        $end = $i;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
        }
    }
    #print "- $#lengths -\n";
    return ($total, \@overlaps, \@lengths);
    @genomeOverlap = (); @overlaps = (); @lengths = (); $total = 0;
}


sub getNonCDSoverlaps {
    # INPUT: ($totalSize, $ref_arrayOverlaps, $ref_arrayLength) = &getNonCDSoverlaps("$size_chr{$species}{$seqname}", \@coordinates); 

    my $positions= (); my @genomeOverlap = (); my @overlaps = (); my @lengths = (); 
    my $pair = ""; my $leftpos = 0; my $rightpos = 0; my $size = 0;
    my $i = 0; my $z = 0;

    $size = shift;
    $positions = shift;
    #$genomeOverlap[$size] = undef;
    
    for ($z= 0; $z < ($size+1); $z++) {
	push @genomeOverlap, "0";
    }

    foreach $pair ( @$positions ) {
        chomp $pair;
	$leftpos = ""; $rightpos = "";
        ($leftpos, $rightpos) = split ":", $pair;
        for ( $i = $leftpos; $i <= $rightpos; $i++ ) { $genomeOverlap[$i] = $genomeOverlap[$i] + 1; }
    }
    $positions= (); 
    my $SizecodingOverlap = 0;
    my $Sizenocoding = 0;
    my $total = 0;
    my $start = 0; my $end = 0;
    my $x = 0; my $y = 0;

    for ( $i = 1; $i <= $#genomeOverlap; $i++ ) {
        if ( $genomeOverlap[$i] == 0 && $Sizenocoding == 0 ) {
            $start = $i;
	    $SizecodingOverlap = 0;
            $Sizenocoding++;
        }
	elsif ( $genomeOverlap[$i] > 0 && $SizecodingOverlap == 0 ) {
            if ( $Sizenocoding > 0 && $i < $#genomeOverlap ) {
	        $end = $i - 1;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
            }
            $SizecodingOverlap++;
            $Sizenocoding  = 0;
        }
	elsif ( $Sizenocoding > 0 && $i == $#genomeOverlap ) {
	        $end = $i;
	        $overlaps[$y++] = "$start:$end";
		$lengths[$x++] = $end-$start;
		$total+= (($end-$start) + 1);
                #print "$start\t$end\n";
                $start = 0;
        }
    }
    #print "- $#lengths -\n";
    return ($total, \@overlaps, \@lengths);
    @genomeOverlap = (); @overlaps = (); @lengths = (); $size = ""; $total = 0;
}




#####################################################################
#                                                                   #
#   CALCULATING SEVERAL STATISTIC POPULATION DESCRIPTORS            #
#                                                                   #
#    - getQuartiles            ---> Quartiles                       #
#    - calculateStdDev         ---> Standard deviation              #
#    - getTypeSizeDistribution ---> Lenght distribution by ranges   #
#                                                                   # 
#####################################################################


sub getQuartiles{
    # INPUT: ($Q1, $Q3, $Q_top, $Q_down, $IQR, $lower_fence, $upper_fence, $median) = &getQuartiles(\@statistics);

    my @values = (); my @Q1 = (); my @Q2 = (); my @Q3 = (); my @percentile_top = (); my @percentile_down = ();
    my ($q_top, $q_down, $q1, $q3, $iqr, $lfence, $ufence, $total, $middle, $medianOne, $medianTwo) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    my ($position, $first_quartile, $third_quartile, $percentile_low, $percentile_high) = (0, 0, 0, 0, 0);

    @values = @{$_[0]}; $total = scalar @values;
    @Q1 = @values; @Q2 = @values; @Q3 = @values; @percentile_top = @values; @percentile_down = @values;

    $first_quartile = (int(0.25*($total + 1))) - 1;
    $third_quartile = (int(0.75*($total + 1))) - 1;
    $percentile_low = (int(0.05*($total + 1))) - 1;
    $percentile_high = (int(0.95*($total + 1))) - 1;
    $q1 = splice @Q1, $first_quartile, 1;
    $q3 = splice @Q3, $third_quartile, 1;
    $q_top = splice @percentile_top, $percentile_high, 1;
    $q_down = splice @percentile_down, $percentile_low, 1;
    $iqr = $q3 - $q1;			# IQR = Interquartile Range
    $lfence = $q1 - (1.5*$iqr);
    $ufence = $q3 + (1.5*$iqr);

    if ( $total % 2 == 0 ) {
	$position = ($total - 2) / 2;
        ($medianOne, $medianTwo) = splice @Q2, $position, 2;
        $middle = int(($medianOne + $medianTwo) / 2);
    }
    else {
        $position = ($total - 1) / 2;
	$middle = splice @Q2, $position, 1;
    }
    #print "$total, $q1/$first_quartile, $q3/$third_quartile, $q_top/$percentile_high, $q_down/$percentile_low, $iqr, $lfence, $ufence, $middle\n";

    return ($q1, $q3, $q_top, $q_down, $iqr, $lfence, $ufence, $middle);

    @values = (); @Q1 = (); @Q2 = (); @Q3 = (); @percentile_top = (); @percentile_down = ();
    ($q_top, $q_down, $q1, $q3, $iqr, $lfence, $ufence, $total, $middle, $medianOne, $medianTwo) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    ($position, $first_quartile, $third_quartile, $percentile_low, $percentile_high) = (0, 0, 0, 0, 0);
}


sub calculateStdDev{
    # INPUT: ($variance, $stddev, $var_filter, $sd_filter, $avg_filter, $null) = &calculateStdDev($AvgCDsGenome, $lower_fence, $upper_fence, \@cds);
    # The Population Standard Deviation: stdDev = sqrt(((sample1 - average)^2 + ... + (sampleN - average)^2)/N)
    # The Sample Standard Deviation: stdDev = sqrt(((sample1 - average)^2 + ... + (sampleN - average)^2)/N -1)

    my ($sumOfSquares, $sumOfSquares_filter, $avg, $average_filter, $sqrt_avg, $sqrt_filter, $variance_avg, $stddev_avg, $variance_filter, $stddev_filter, $outlier_down, $outlier_top) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    my ($x, $y, $z, $n_all, $n_filter, $sum_all, $sum_filter, $sumsquared_diff, $sumsquared_diff_filter, $avg_all) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    my $values = (); my @values_filter = (); my %hash_sizescount = ();

    $avg = shift; $outlier_down = shift; $outlier_top = shift;
    $values = shift;
    #print "@$values\n*************************\n";

    foreach $x ( @$values ) {
        chomp $x;
	$x =~ s/\s|\t//gis;
	$hash_sizescount{all}{$x}++;
	$sumsquared_diff+= (($x - $avg) ** 2);
        $sum_all+= $x;
        $n_all++;
	$sumOfSquares+= $x * $x;

	if ( $x >= $outlier_down && $x <= $outlier_top ) {
	    $sum_filter += $x;
	    $n_filter++;
	    $values_filter[$z++] = $x;
	    $hash_sizescount{filter}{$x}++;
	}
    }
    if ( $n_filter > 1 ) {
        $average_filter = $sum_filter/$n_filter;

        foreach $y ( @values_filter ) {
    	    $sumsquared_diff_filter+= (($y - $average_filter) ** 2);
	    $sumOfSquares_filter+= $y * $y;
        }
    }

    if ( $n_all > 1 && $n_filter > 1 ) {
        $avg_all = sprintf("%.2f", ($sum_all / $n_all));
        $sqrt_avg = sqrt($sumsquared_diff / $n_all);
	$sqrt_filter = sqrt($sumsquared_diff_filter / ($n_filter - 1));
        $variance_avg = sprintf("%.2f", (($sumOfSquares - (($sum_all * $sum_all) / $n_all)) / $n_all));
        $stddev_avg = sprintf("%.2f", (sqrt($variance_avg)));
	$variance_filter = sprintf("%.2f", (($sumOfSquares_filter - (($sum_filter * $sum_filter) / $n_filter)) / ($n_filter - 1)));
        $stddev_filter = sprintf("%.2f", (sqrt($variance_filter)));
	$average_filter = sprintf("%.2f", ($average_filter));
    }
    else { $variance_avg = 0; $stddev_avg = 0; $variance_filter = 0; $stddev_filter = 0; $average_filter = 0; %hash_sizescount = (); $avg_all = 0; }
    #foreach (keys %hash_sizescount) { print "$_ = $hash_sizescount{$_}\n"; }
    #print "new standard deviation : $sqrt_avg = $sqrt_filter \t old standard deviation : $stddev_avg = $stddev_filter\n";
    #print "$variance_avg => $stddev_avg => $variance_filter => $stddev_filter => $average_filter => $avg_all\n";

    return ($variance_avg, $stddev_avg, $variance_filter, $stddev_filter, $average_filter, \%hash_sizescount, $avg_all);

    $values = (); @values_filter = (); %hash_sizescount = ();
    ($sumOfSquares, $sumOfSquares_filter, $avg, $average_filter, $sqrt_avg, $sqrt_filter, $variance_avg, $stddev_avg, $variance_filter, $stddev_filter, $outlier_down, $outlier_top) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}


sub getTypeSizeDistribution{
    # INPUT: ($ref_hashDistribution, $ref_hashRangesATGC, $summatory) = &getTypeSizeDistribution(\@atgcs);

    my $values = (); my @ranges = (); my %hashByRanges = (); my %hashByRangesATGC = ();
    my $total = 0; my $data = "";
    my ($start, $end) = (0, 0);
    my ($nucleotides, $range, $x, $at, $gc, $n) = ("", "", 0, 0, 0, 0);

    $values = shift;
    $total = scalar @$values;
    @ranges = (      "1-15",           "16-50",             "51-100",           "101-250",              "251-500",             "501-1000", 
	          "1001-5000",       "5001-10000",       "10001-50000",       "50001-100000",       "100001-500000",       "500001-1000000", 
	       "1000001-5000000", "5000001-10000000", "10000001-50000000", "50000001-100000000", "100000001-500000000", "500000001-1000000000");

    foreach $range ( @ranges ) {
        chomp $range;
	($start, $end) = (0, 0);

	($start, $end) = split "-", $range;
	$hashByRangesATGC{$end}{AT} = 0; $hashByRangesATGC{$end}{GC} = 0; $hashByRangesATGC{$end}{N} = 0;
	$hashByRanges{$end} = 0;

        foreach $data ( @$values ) {
            chomp $data;							## "$lengthExon=$atExon:$gcExon:$nExon";
	    ($nucleotides, $range, $x, $at, $gc, $n) = ("", "", 0, 0, 0, 0);

	    ($x, $nucleotides) = split "=", $data;
	    ($at, $gc, $n) = split ":", $nucleotides;
	    $at =~ s/\s|\t//gis; $gc =~ s/\s|\t//gis; $n =~ s/\s|\t//gis;
	    #print "$x ---> $nucleotides = $start <=> $end <=> $at, $gc, $n\n";

            if ( $x >= $start && $x <= $end ) { 
	        $hashByRanges{$end}++;
		$hashByRangesATGC{$end}{AT}+= $at;
		$hashByRangesATGC{$end}{GC}+= $gc;
		$hashByRangesATGC{$end}{N}+= $n;
	    }
	}
    }
    return (\%hashByRanges, \%hashByRangesATGC, $total);
    $values = (); @ranges = (); %hashByRanges = (); %hashByRangesATGC = ();
}




#################################################################
#                                                               #
#   ESTIMATING THE QUALITY OF PROTEIN-CODING GENE ANNOTATIONS   #
#   Intron length modulo 3, according to Roy and Penny (2007)   #
#                                                               # 
#################################################################


sub getIntronModulo3Distributions{
    #INPUT &getIntronModulo3Distributions("$genera $specie", $outfiles, $single, 
    #\%NumberBySpecies, \%NumberByChromosomes, \%NumberByStrands, 
    #\%IntronClassBySpecies, \%IntronClassByChromosomes, \%IntronClassByStrands);

    my $sps = ""; my $dir = ""; my $process = ""; my $plusORminus = "";
    my %number_sp = (); my %number_chr = (); my %number_str = (); my %class_sp = (); my %class_chr = (); my %class_str = ();
    my $total_3n = 0; my $total_3n1 = 0; my $total_3n2 = 0;
    my $fraction_3n = 0; my $fraction_3n1 = 0; my $fraction_3n2 = 0;
    my $excess_3n = 0; my $excess_3n12 = 0;

    $sps = shift; $dir = shift; $process = shift;
    %number_sp = %{$_[0]}; %number_chr = %{$_[1]}; %number_str = %{$_[2]}; %class_sp = %{$_[3]}; %class_chr = %{$_[4]}; %class_str = %{$_[5]};

    if ( $process =~ /[NnFf]/ ) {
        open STATSTWO_STRAND, ">>$dir/total.intron_types.strands.table.species.txt" || die "Cannot open $dir/total.intron_types.strands.table.species.txt\n";
        open STATSTWO_GENOME, ">>$dir/total.intron_types.table.species.txt" || die "Cannot open total.intron_types.table.species.txt\n";
    }
    if ( $process =~ /[YyTt]/ ) {
        open STATSTWO_STRAND, ">>$dir/total.$sps.intron_types.strands.table.txt" || die "Cannot open $dir/total.$sps.intron_types.strands.table.txt\n";
        open STATSTWO_GENOME, ">>$dir/total.$sps.intron_types.table.txt" || die "Cannot open $dir/total.$sps.intron_types.table.txt\n";
    }

    if ( $number_sp{$sps}{introns} ) {
        $total_3n = ""; $total_3n1 = ""; $total_3n2 = ""; $excess_3n = ""; $excess_3n12 = "";

        $total_3n = sprintf("%.3f",(($class_sp{$sps}{'3n'} / $number_sp{$sps}{introns})));
        $total_3n1 = sprintf("%.3f",(($class_sp{$sps}{'3n+1'} / $number_sp{$sps}{introns})));
        $total_3n2 = sprintf("%.3f",(($class_sp{$sps}{'3n+2'} / $number_sp{$sps}{introns})));
        $excess_3n = sprintf("%.3f",($total_3n - (($total_3n1 + $total_3n2) / 2)));
        $excess_3n12 = sprintf("%.3f",($total_3n1 - $total_3n2));
        print STATSTWO_GENOME "$sps\t$number_sp{$sps}{introns}\t$total_3n\t$total_3n1\t$total_3n2\t$excess_3n\t$excess_3n12\n";
        #print "$sps\t$number_sp{$sps}{introns}\t$total_3n\t$total_3n1\t$total_3n2\t$excess_3n\t$excess_3n12\n";
    }
    else { print STATSTWO_GENOME "$sps\t0\t\t0\t0\t0\t0\n"; }

    for $plusORminus ( sort keys %{ $class_str{$sps} } ) {
        if ( $number_str{$sps}{introns}{$plusORminus} ) {
            $total_3n = ""; $total_3n1 = ""; $total_3n2 = ""; $excess_3n = ""; $excess_3n12 = "";

            $fraction_3n = sprintf("%.3f",($class_str{$sps}{$plusORminus}{'3n'} / $number_str{$sps}{introns}{$plusORminus}));
	    $fraction_3n1 = sprintf("%.3f",($class_str{$sps}{$plusORminus}{'3n+1'} / $number_str{$sps}{introns}{$plusORminus}));
	    $fraction_3n2 = sprintf("%.3f",($class_str{$sps}{$plusORminus}{'3n+2'} / $number_str{$sps}{introns}{$plusORminus}));
	    $excess_3n = sprintf("%.3f",($fraction_3n - (($fraction_3n1 + $fraction_3n2) / 2)));
	    $excess_3n12 = sprintf("%.3f",($fraction_3n1 - $fraction_3n2));
	    print STATSTWO_STRAND "$sps\t$plusORminus\t$number_str{$sps}{introns}{$plusORminus}\t$fraction_3n\t$fraction_3n1\t$fraction_3n2\t$excess_3n\t$excess_3n12\n";
	    #print "$sps\t$plusORminus\t$number_str{$sps}{introns}{$plusORminus}\t$fraction_3n\t$fraction_3n1\t$fraction_3n2\t$excess_3n\t$excess_3n12\n"; 
	}
	else { print STATSTWO_STRAND "$sps\t$plusORminus\t0\t0\t0\t0\t0\t0\n"; }
    }
    $sps = ""; $dir = ""; $process = "";  $sps = "";
    %number_sp = (); %number_chr = (); %number_str = (); %class_sp = (); %class_chr = (); %class_str = ();
}




#################################################################
#                                                               #
#   PREPARING THE FORMAT OF COORDINATES FOR OUTPUT FILES        #
#   Include quartile estimations to coodinates                  #
#                                                               # 
#################################################################


sub addQuartilesToCoordinates {
    # INPUT: &addQuartilesToCoordinates($outfile_gff, $type, $Q1, $Q3, $Q_top, $Q_down, $lower_fence, $upper_fence, $median, \%all_outfile);

    my %coordinatesOut = ();
    my $outall = ""; my $class = ""; my $i = ""; my $fastaname = ""; my $all = ""; my $length = 0; my $plusORminus = "";
    my $Qone = 0; my $Qthree = 0; my $Qtop = 0; my $Qdown = 0; my $LOWoutlier = 0; my $UPoutlier = 0; my $Qtwo = 0;

    $outall = shift; $class = shift; $Qone = shift; $Qthree = shift; $Qtop = shift; $Qdown = shift; $LOWoutlier = shift; $UPoutlier = shift; $Qtwo = shift;
    %coordinatesOut = %{$_[0]};

    open ALL, ">>$outall" || die "Cannot open $outall\n";

    for $fastaname ( sort keys %coordinatesOut ) {
	for $plusORminus ( sort keys %{ $coordinatesOut{$fastaname}{$class} } ) {
            for $i ( 0 .. $#{ $coordinatesOut{$fastaname}{$class}{$plusORminus} } ) {
		#print "$coordinatesOut{$fastaname}{$class}{$plusORminus}[$i]\t=\t$Qone, $Qtwo, $Qthree, $Qtop, $Qdown, $LOWoutlier, $UPoutlier\n";
	        ($all, $length) = split "=", $coordinatesOut{$fastaname}{$class}{$plusORminus}[$i];
		if ( $length < $Qone ) { print ALL "$all\t$length\tQ0:<$Qone\n"; }
		if ( $length >= $Qone && $length < $Qtwo ) { print ALL "$all\t$length\tQ1-Q2:$Qone-$Qtwo\n"; }
		if ( $length > $Qtwo && $length <= $Qthree ) { print ALL "$all\t$length\tQ2-Q3:$Qtwo-$Qthree\n"; }
		if ( $length > $Qthree ) { print ALL "$all\t$length\tQ4:>$Qthree\n"; }
	    }
	}
    }
    %coordinatesOut = (); $outall = ""; $class = ""; $Qone = 0; $Qthree = 0; $Qtop = 0; $Qdown = 0; $LOWoutlier = 0; $UPoutlier = 0; $Qtwo = 0;
}




#################################################################
#                                                               #
#   OPENING THE OUTPUT FILES IN THE CORRESPONDING DIRECTORIES   #
#                                                               # 
#################################################################


sub openOutfilesSeveralOption {
    # INPUT: if ( $single=~ /[NnFf]/ ) { &openOutfilesSeveralOption($outfiles); }

    my $dir = ""; $dir = shift;

    open SPECIES, ">>$dir/total.summary.genome.txt" || die "Cannot open $dir/total.summary.genome.txt\n";
    print SPECIES "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tCDS_AVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\n";

    open STATSEXONS_CHROMOSOME, ">>$dir/total.exons.chromosomes.table.txt" || die "Cannot open $dir/total.exons.chromosomes.table.txt\n";
    open STATSEXONS_STRANDS, ">>$dir/total.exons.strands.table.txt" || die "Cannot open $dir/total.exons.strands.table.txt\n";
    open STATSEXONS_GENOME, ">>$dir/total.exons.table.txt" || die "Cannot open $dir/total.exons.table.txt\n";

    open STATSINTRONS_CHROMOSOME, ">>$dir/total.introns.chromosomes.table.txt" || die "Cannot open $dir/total.introns.chromosomes.table.txt\n";
    open STATSINTRONS_STRANDS, ">>$dir/total.introns.strands.table.txt" || die "Cannot open $dir/total.introns.strands.table.txt\n";
    open STATSINTRONS_GENOME, ">>$dir/total.introns.table.txt" || die "Cannot open $dir/total.introns.table.txt\n";

    open STATSINTERGENICS_CHROMOSOME, ">>$dir/total.intergenics.chromosomes.table.txt" || die "Cannot open $dir/total.intergenics.chromosomes.table.txt\n";
    open STATSINTERGENICS_STRANDS, ">>$dir/total.intergenics.strands.table.txt" || die "Cannot open $dir/total.intergenics.strands.table.txt\n";
    open STATSINTERGENICS_GENOME, ">>$dir/total.intergenics.table.txt" || die "Cannot open $dir/total.intergenics.table.txt\n";

    print STATSEXONS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSEXONS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSEXONS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tTOTAL_NUMBER\tSTANDARD_DENSITY\tDENSITY_MIN\tDENSITY_MAX\tDENSITY_VAR\tDENSITY_SD\tDENSITY_EBAR_BEGIN\tDENSITY_EBAR_END\tDENSITY_FILTER\tDENSITY_FILTER_VAR\tDENSITY_FILTER_SD\tDENSITY_FILTER_EBARBEGIN\tDENSITY_FILTER_EBAREND\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";

    print STATSINTRONS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tCDS%\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSINTRONS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tCDS%\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVE\tLOG-MEAN\n";
    print STATSINTRONS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_with_INTRONS%\tTOTAL_NUMBER\tEXON_DEFINITION%\tEXONFILTER_DEFINITION%\tSTANDARD_DENSITY\tWEIGHTED_DENSITY\tRATIO_INTRONS:EXONS\tDENSITY_MIN\tDENSITY_MAX\tDENSITY_VAR\tDENSITY_SD\tDENSITY_EBAR_BEGIN\tDENSITY_EBAR_END\tDENSITY_FILTER\tDENSITY_FILTER_VAR\tDENSITY_FILTER_SD\tDENSITY_FILTER_EBARBEGIN\tDENSITY_FILTER_EBAREND\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";

    print STATSINTERGENICS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tGENOME_CDS_%\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";
    print STATSINTERGENICS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";
    print STATSINTERGENICS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tGENOME_CDS_%\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";

    open STATSTWO_STRAND, ">>$dir/total.intron_types.strands.table.txt" || die "Cannot open $dir/total.intron_types.strands.table.txt\n";
    open STATSTWO_GENOME, ">>$dir/total.intron_types.table.txt" || die "Cannot open $dir/total.intron_types.table.txt\n";

    print STATSTWO_STRAND "SPECIES\tSTRAND\tINTRONS\t3n\t3n+1\t3n+2\tEXCESS 3n\tEXCESS 3n+1\n";
    print STATSTWO_GENOME "SPECIES\tINTRONS\t3n\t3n+1\t3n+2\tEXCESS 3n\tEXCESS 3n+1\n";

    open EXONSRANGE_GENOME, ">>$dir/total.exons_ranges.table.txt" || die "Cannot open $dir/total.exons_ranges.table.txt\n";
    open INTRONSRANGE_GENOME, ">>$dir/total.introns_ranges.table.txt" || die "Cannot open $dir/total.introns_ranges.table.txt\n";
    open INTERGENICSRANGE_GENOME, ">>$dir/total.intergenics_ranges.table.txt" || die "Cannot open $dir/total.intergenics_ranges.table.txt\n";

    print EXONSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";
    print INTRONSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";
    print INTERGENICSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";

    open EXONSATGC_GENOME, ">>$dir/total.exons_ranges_atgc.table.txt" || die "Cannot open $dir/total.exons_ranges_atgc.table.txt\n";
    open INTRONSATGC_GENOME, ">>$dir/total.introns_ranges_atgc.table.txt" || die "Cannot open $dir/total.introns_ranges_atgc.table.txt\n";
    open INTERGENICSATGC_GENOME, ">>$dir/total.intergenics_ranges_atgc.table.txt" || die "Cannot open $dir/total.intergenics_ranges_atgc.table.txt\n";

    print EXONSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
    print INTRONSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
    print INTERGENICSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
}


sub openOutfilesSingleOption {
    # INPUT: if ( $single=~ /[YyTt]/ ) { &openOutfilesSingleOption($outfiles, $species); }

    my $dir = ""; my $sps = "";
    $dir = shift; $sps = shift;

    open SPECIES, ">>$dir/total.$sps.summary.genome.txt" || die "Cannot open $dir/total.$sps.summary.genome.txt\n";
    print SPECIES "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tCDS_AVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\n";

    open STATSEXONS_CHROMOSOME, ">>$dir/total.$sps.exons.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.exons.chromosomes.table.txt\n";
    open STATSEXONS_STRANDS, ">>$dir/total.$sps.exons.strands.table.txt" || die "Cannot open $dir/total.$sps.exons.strands.table.txt\n";
    open STATSEXONS_GENOME, ">>$dir/total.$sps.exons.table.txt" || die "Cannot open $dir/total.$sps.exons.table.txt\n";

    open STATSINTRONS_CHROMOSOME, ">>$dir/total.$sps.introns.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.introns.chromosomes.table.txt\n";
    open STATSINTRONS_STRANDS, ">>$dir/total.$sps.introns.strands.table.txt" || die "Cannot open $dir/total.$sps.introns.strands.table.txt\n";
    open STATSINTRONS_GENOME, ">>$dir/total.$sps.introns.table.txt" || die "Cannot open $dir/total.$sps.introns.table.txt\n";

    open STATSINTERGENICS_CHROMOSOME, ">>$dir/total.$sps.intergenics.chromosomes.table.txt" || die "Cannot open $dir/total.$sps.intergenics.chromosomes.table.txt\n";
    open STATSINTERGENICS_STRANDS, ">>$dir/total.$sps.intergenics.strands.table.txt" || die "Cannot open $dir/total.$sps.intergenics.strands.table.txt\n";
    open STATSINTERGENICS_GENOME, ">>$dir/total.$sps.intergenics.table.txt" || die "Cannot open $dir/total.$sps.intergenics.table.txt\n";

    print STATSEXONS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSEXONS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSEXONS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tTOTAL_NUMBER\tSTANDARD_DENSITY\tDENSITY_MIN\tDENSITY_MAX\tDENSITY_VAR\tDENSITY_SD\tDENSITY_EBAR_BEGIN\tDENSITY_EBAR_END\tDENSITY_FILTER\tDENSITY_FILTER_VAR\tDENSITY_FILTER_SD\tDENSITY_FILTER_EBARBEGIN\tDENSITY_FILTER_EBAREND\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";

    print STATSINTRONS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tCDS%\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";
    print STATSINTRONS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tCDS%\tSTANDARD_DENSITY\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVE\tLOG-MEAN\n";
    print STATSINTRONS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_with_INTRONS%\tTOTAL_NUMBER\tEXON_DEFINITION%\tEXONFILTER_DEFINITION%\tSTANDARD_DENSITY\tWEIGHTED_DENSITY\tRATIO_INTRONS:EXONS\tDENSITY_MIN\tDENSITY_MAX\tDENSITY_VAR\tDENSITY_SD\tDENSITY_EBAR_BEGIN\tDENSITY_EBAR_END\tDENSITY_FILTER\tDENSITY_FILTER_VAR\tDENSITY_FILTER_SD\tDENSITY_FILTER_EBARBEGIN\tDENSITY_FILTER_EBAREND\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tWEIGHTAVG_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tWEIGHTAVG_VAR\tWEIGHTAVG_SD\tWEIGHTAVG_EBARBEGIN\tWEIGHTAVG_EBAREND\tWEIGHTAVG_VARFILTER\tWEIGHTAVG_SDFILTER\tWEIGHTAVG_EBARBEGIN_FILTER\tWEIGHTAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG_AVG_SIZE\tLOG_WAVG_SIZE\n";

    print STATSINTERGENICS_CHROMOSOME "SPECIES\tCHROMOSOME\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tGENOME_CDS_%\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";
    print STATSINTERGENICS_STRANDS "SPECIES\tGENOME_MBS\tSTRAND\tCDS_NUMBER\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tCONTENT_GENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tAVG_VAR\tAVG_SD\tAVG_EBARBEGIN\tAVG_EBAREND\tAVG_VARFILTER\tAVG_SDFILTER\tAVG_EBARBEGIN_FILTER\tAVG_EBAREND_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";
    print STATSINTERGENICS_GENOME "SPECIES\tGENOME_MBS\tCDS_NUMBER\tCDS_CONTENT_MBS\tGENOME_CDS_%\tTOTAL_NUMBER\tFULL_LENGTH_MBS\tFULL_LENGTH_NTS\tCONTENT_MBS\tCONTENT_NTS\tGENOME_%\tMIN_SIZE\tMAX_SIZE\tAVERAGE_SIZE\tFILTERAVG_SIZE\tNORMALISEDAVG_SIZE\tSIZE_Q_DOWN\tSIZE_Q1\tSIZE_Q2-MEDIAN\tSIZE_Q3\tSIZE_Q_UP\tSIZE_LOWER_FENCE\tSIZE_UPPER_FENCE\tVARIENCE\tSD\tERRORBAR_BEGIN\tERRORBAR_END\tVAR_FILTER\tSD_FILTER\tEBAR_BEGIN_FILTER\tEBAR_END_FILTER\tFREQ_A%\tFREQ_T%\tFREQ_G%\tFREQ_C%\tFREQ_Ns%\tFREQ_AT%\tFREQ_GC%\tLOG_GENOME_MBS\tLOG_GENOME_CONTENT\tLOG-AVG\n";

    open STATSTWO_STRAND, ">>$dir/total.$sps.intron_types.strands.table.txt" || die "Cannot open $dir/total.$sps.intron_types.strands.table.txt\n";
    open STATSTWO_GENOME, ">>$dir/total.$sps.intron_types.table.txt" || die "Cannot open $dir/total.$sps.intron_types.table.txt\n";

    print STATSTWO_STRAND "SPECIES\tSTRAND\tINTRONS\t3n\t3n+1\t3n+2\tEXCESS 3n\tEXCESS 3n+1\n";
    print STATSTWO_GENOME "SPECIES\tINTRONS\t3n\t3n+1\t3n+2\tEXCESS 3n\tEXCESS 3n+1\n";

    open EXONSRANGE_GENOME, ">>$dir/total.$sps.exons_ranges.table.txt" || die "Cannot open $dir/total.$sps.exons_ranges.table.txt\n";
    open INTRONSRANGE_GENOME, ">>$dir/total.$sps.introns_ranges.table.txt" || die "Cannot open $dir/total.$sps.introns_ranges.table.txt\n";
    open INTERGENICSRANGE_GENOME, ">>$dir/total.$sps.intergenics_ranges.table.txt" || die "Cannot open $dir/total.$sps.intergenics_ranges.table.txt\n";

    print EXONSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";
    print INTRONSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";
    print INTERGENICSRANGE_GENOME "SPECIES (number)\tGENOME_MBS\t15\t50\t100\t250\t500\t1000\t5000\t10000\t50000\t100000\t500000\t1000000\t5000000\t10000000\t50000000\t100000000\t500000000\t1000000000\n";

    open EXONSATGC_GENOME, ">>$dir/total.$sps.exons_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.exons_ranges_atgc.table.txt\n";
    open INTRONSATGC_GENOME, ">>$dir/total.$sps.introns_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.introns_ranges_atgc.table.txt\n";
    open INTERGENICSATGC_GENOME, ">>$dir/total.$sps.intergenics_ranges_atgc.table.txt" || die "Cannot open $dir/total.$sps.intergenics_ranges_atgc.table.txt\n";

    print EXONSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
    print INTRONSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
    print INTERGENICSATGC_GENOME "SIZE_NTS_RANGE\tNUMBER\tFREQ_GC%\tFREQ_AT%\tFREQ_N%\n";
}
