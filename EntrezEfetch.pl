#! /bin/env perl
# download sequences based on the results of Entrez or the list of sequences IDs
# Updated: 2014-09-09
# The usage information is available by specifying '-h' or '--help'
# Updated: 2014-08-13

# Written by Sishuo Wang from the Department of Botany, the University of British Columbia
# Please write e-mail to sishuowang@hotmail.ca if you have any question and/or suggestion. Your help will be highly appreciated!

############################################################################################
use File::Basename;
use Bio::DB::EUtilities;
use Bio::DB::Query::GenBank;
use Getopt::Long;
use strict;

my ($db, $searching_word, $retmax, $outfile, $outdir, $format, $ids_list_file, $is_combine, $num_of_procs, $is_clean, $is_no_blank_line);
my (@ids);
my ($gb, $seqio, $query, $is_entrez, $force, $help);

$format='fasta';
$retmax=100;
$num_of_procs=1;
$is_no_blank_line=0;

GetOptions(
    'searching_word|key_word=s'	=>  \$searching_word,
    'db=s'		    =>	\$db,
    'retmax=s'		=>	\$retmax,
    'outfile=s'		=>	\$outfile,
    'outdir=s'		=>	\$outdir,
    'format=s'		=>	\$format,
    'id|gi=s'		=>	\@ids,
    'ids_list=s'	=>	\$ids_list_file,
    'combine!'		=>	\$is_combine,
    'procs=s'		=>	\$num_of_procs,
    'clean!'		=>	\$is_clean,
    'force!'		=>	\$force,
    'h|help!'		=>	\$help,
    'no_blank_line!'    =>      \$is_no_blank_line,
) || die "illegal params!";


&USAGE() if $help;
warn "Neither outfile nor outdir was specified!\n" and &USAGE() if not $outfile and not $outdir;
warn "num_of_procs should be an integer between 1 and 4" and &USAGE() if $num_of_procs !~ /^[1-4]$/;
if ($outdir){
    -d $outdir and `rm -rf $outdir` if $force;
    `mkdir -p $outdir`;
}
$is_entrez=1 if $searching_word;


###########################################################
if ($is_entrez){
    @ids=@{&entrez($db, $searching_word, $retmax)};
}
else{
    @ids=@{&read_ids_list($ids_list_file)} if $ids_list_file;
}

print "Downloading ......\n";
my ($lower_limit, $upper_limit);
$upper_limit=-1;
while($upper_limit<$#ids){
    $lower_limit=$upper_limit+1;
    $upper_limit=($lower_limit+499 >= $#ids ? $#ids : $lower_limit+499);
    my @new_ids=@ids[$lower_limit..$upper_limit];
    print eval{$lower_limit+1}."\t".($upper_limit+1)."\n";

    my $num_per_Fen = int(scalar(@new_ids)/$num_of_procs);
    my ($start_k, $end_k, $i);
    my (@children_pids);
    $end_k=-1;
    for ($i=1; $i<=$num_of_procs; $i++){
        $start_k=$end_k+1;
        $end_k=($start_k+$num_per_Fen >= $#new_ids ? $#new_ids : $start_k+$num_per_Fen);
        my $pid = fork();
        if ($pid == 0){
            $outfile=$outdir."/"."outfile.".join("-",($lower_limit+1,$upper_limit+1)).".proc$i";
            &download([@new_ids[$start_k..$end_k]], $outfile);
            my $basename=basename($outfile);
            #print $outfile."\n";
            #&download([@new_ids], $outfile);
            exit 0;
        }
        else{
            push @children_pids, $pid;
        }
    }
    foreach my $child (@children_pids) {
        waitpid($child, 0); 
    }
}


if (not $outfile){
    $outfile = $outdir."/"."final_output";
}
&after_analyses($outfile);

###########################################################
sub entrez{
    my ($db, $searching_word, $retmax) = @_;
    my $factory = Bio::DB::EUtilities->new(-eutil  => 'esearch',
                                           -db     => $db,
                                           -term   => $searching_word,
                                           -email  => 'mymail@foo.bar',
                                           -retmax => $retmax);

    my @ids = $factory->get_ids;
    return (\@ids);
}

sub read_ids_list{
    my @ids;
    my $ids_list_file=shift;
    open(my $IN, '<', $ids_list_file) || die "failed to open ids_list_file $ids_list_file";
    while(<$IN>){
        chomp;
        push @ids, $_;
    }
    return(\@ids);
}


sub download{
    my ($new_ids_aref, $outfile) = @_;
    #do {print $_."\n"} foreach @$new_ids_aref;
    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                           -db      => $db,
                                           -rettype => $format,
                                           -id      => $new_ids_aref);
    $factory->get_Response(-file => $outfile);
    return();
    my $guesser = new Bio::Tools::GuessSeqFormat( -file => $outfile);
    my $format  = $guesser->guess;
    my $seqin = Bio::SeqIO->new(-file   => $outfile,
                                -format => $format,
                                );

    while (my $seq = $seqin->next_seq) {
        print $seq->display_id."\n";
    }
}

sub after_analyses{
    no strict;
    my $combine_outfile = shift;
    $combine_outfile="$outdir/combined_seqs.fas" if not $combine_outfile;
    if ($is_combine and $outdir){
        `for i in $outdir/*proc*; do sed -i '/^$/d' $i; done` if $is_no_blank_line;
        `cat $outdir/*proc* > $combine_outfile`;
        if ($is_clean){
	        foreach(glob("$outdir/*proc*")){
		        `rm $_`;
	        }
	    }
    }
}

sub USAGE{
    my $basename=basename $0;
    print "The usage of $0:\n";
    print "Note: Bioperl needs to be installed in advance\n";
    print "perl $basename [Options]\n";
    print<<EOF
    --key_word	key_word used in Entrez
    --db	database to search, e.g. nucleotide, nucest, protein, etc.
    --retmax	maximun number of sequences to retrieve
		default: 100
    --outfile	output file (only functional when --combine is specified)
    --outdir	outdir
    --format	format of sequences to download
		default: fasta
    --id|--gi	id to search (GI)
    --ids_list	file containing the list of ids (GI) of sequences to download
    --combine	combine results into one file (no argument required)
		default: NO
    --proc	number of processes (1-4)
		default: 1
    --clean	clean separated files (no argument required)
		default: NO
    --force	force to remove outdir if it exists
		default: NO
    -h|--help	see usage
EOF
;
    print "If you have any question or suggestion, please write e-mail to sishuowang\@hotmail.ca. Thanks!\n\n";
    exit(1);
}

