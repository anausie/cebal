use strict;
use warnings;
use File::Path;
############
#Main start#
############

########## PARAMS #########
# - Input File
# - Data type
###########################################################
# Input File: 5 columns separated by TABs.
# - 1st Column: Path to the SFF file
# - 2nd Column: Key Tag sequence
# - 3rd Column: Adaptor sequence
# - 4th Column: 0 if key tag is present; 1, otherwise
# - 5th Column: 0 if adaptor is present; 1, otherwise.
###########################################################
# Data type:
# - 0 for transcriptomic data; 1 for genomic data
###########################################################

if ( @ARGV != 1 )
{
  print "Error arguments!\n";
  print "You must an input file.";
  print "Input File: 5 columns separated by TABs.";
  print "- 1st Column: Path to the SFF file";
  print "- 2nd Column: Key Tag sequence";
  print "- 3rd Column: Adaptor sequence";
  print " - 4th Column: 0 if key tag is present; 1, otherwise";
  print " - 5th Column: 0 if adaptor is present; 1, otherwise.";
  exit(-1)
}

my $listFiles=$ARGV[0]; # list of SFF files with key tag and adapter definitions.
#my $dataType=$ARGV[1];


my %files=(); # hash with the SFF file name and the length of the keyTag+adaptor that will be the start trim point.
my $initIndex=0;
my $fasta="";
my $finalTP="";
# Reading the input file to stablish the starting trim point of each SFF file (all the reads in the same files has the same value).
open (IH, "<", $listFiles) or die "Could not open file '$listFiles' $!";
while (my $line = <IH>){
  chomp($line);
  my @array = split(/\s+/, $line);
  my $name=$array[0];
  my $tag_seq=$array[1];
  my $adaptor_seq=$array[2];
  my $hasTag=$array[3]; # 0-> YES 1-> NO
  my $hasAdaptor=$array[4]; # 0-> YES 1-> NO
  if($hasAdaptor == 0 && $hasTag==0){ #0-> YES 1-> NO
    $initIndex=length($tag_seq)+length($adaptor_seq);
        }
  elsif($hasAdaptor == 1 && $hasTag==0){
    $initIndex=length($adaptor_seq);
        }
  elsif($hasAdaptor == 0 && $hasTag==1){
    $initIndex=length($tag_seq);
        }
   $files{$name} = $initIndex+1;      # hash, using variables
  }
close(IH);
## Create the working directory
my $workingDirectory="preProcessingTmp";
mkdir $workingDirectory;

## For each SFF file do the follow steps
# 1st: Remove the adapters.
# 2nd: Convert to fasta
# 3rd: Execute SeqClean
# 4rt: Trim the reads within the SFF file using new trim points
# 5th: Convert to fasta & qual


for my $key ( keys %files ) {
  my $starting = $files{$key};
  my $ending=0; ## sfffile: A value of 0 specifies that end of the read should be used
  ## 1st STEP: Removing adapters
  my $ind = rindex($key, '/');
  my $nameFile=substr($key, $ind+1);
  print $nameFile."\n";
  my $adaptorTP=$workingDirectory."/".$nameFile.".adapTP";
  open(my $fh, '>', $adaptorTP) or die "Could not open file '$adaptorTP' $!";
  my @list = split $/, `sffinfo -a $key`; # Get the Reads ID to create the TP file
  foreach(@list){
    print $fh "$_\t$starting\t$ending\n";
    }
  close($fh);
  my $outputAdaptorsTP=$workingDirectory."/".$nameFile.".adap.TRIM";
  print "$outputAdaptorsTP\n";
  my $removeAdaptorsSFF="sfffile -o $outputAdaptorsTP -tr $adaptorTP $key";
  print "COMMAND: $removeAdaptorsSFF\n";
  system($removeAdaptorsSFF);
  ## 2nd STEP: Convert 2 fasta
  my $adaptorsFasta=$workingDirectory."/".$nameFile.".fna";
  my $sff2fasta="sffinfo -s $outputAdaptorsTP > $adaptorsFasta";
  print "COMMAND: $sff2fasta \n";
  system($sff2fasta);

  ## 3rd STEP: Execute SeqClean
  my $outputSeqClean=$workingDirectory."/".$nameFile.".cleaned";
  my $seqclean="seqclean $adaptorsFasta -n 10000 -o $outputSeqClean";
  print "COMMAND: ".$seqclean."\n";
  system($seqclean);

  ## 4rt STEP: Trim the reads
  #--- SeqClean Trim Points
  my $outputCln=$nameFile.".fna.cln";
  my $outputValidIDs=$workingDirectory."/".$nameFile.".newIDS";
  my $getIDS="cat $outputCln | grep -vE \"(short|low_qual|dust)\" | cut -f 1 > $outputValidIDs";
  print "COMMAND: $getIDS \n";
  system($getIDS);
  my $outputIDS=$workingDirectory."/".$nameFile.".IDS";
  my $getReadsIDS="sfffile -o $outputIDS -i $outputValidIDs $key";
  print "COMMAND: $getReadsIDS \n";
  system($getReadsIDS);
  my $outputSeqCleanTP=$workingDirectory."/".$nameFile.".seqCleanTP";
  my $getTrimPoints="cat $outputCln | grep -vE \"(short|low_qual|dust)\" | cut -f 1,3,4 > $outputSeqCleanTP";
  print "COMMAND: $getTrimPoints \n";
  system($getTrimPoints);
  #--- Merge Trim Points: AdaptTP + seqCleanTP;
  $finalTP=$workingDirectory."/".$nameFile.".TP";
  open($fh, '>', $finalTP) or die "Could not open file '$finalTP' $!";
  open(IH, '<', $outputSeqCleanTP) or die "Could not open file '$outputSeqCleanTP' $!";
  while (my $line = <IH>){
    chomp($line);
    my @array = split(/\s+/, $line);
    my $readID=$array[0];
    my $ini=$array[1];
    my $end=$array[2];
    $ini=$ini+$starting-1;
    $end=$end+$starting-1;
    print $fh "$readID\t$ini\t$end\n";
    }
  close(IH);
  close($fh);
  #--- Apply the new Trim Points
  my $outputTRIM=$workingDirectory."/".$nameFile.".TRIM";
  my $getReadsTRIM="sfffile -o $outputTRIM -tr $finalTP $outputIDS";
  print "COMMAND: $getReadsTRIM \n";
  system($getReadsTRIM);
  # 5th STEP: Convert to fasta+qual
  $fasta=$workingDirectory."/".$nameFile.".TRIM.fasta";
  my $qual=$workingDirectory."/".$nameFile.".TRIM.qual";
  $sff2fasta="sffinfo -s $outputTRIM > $fasta";
  print "COMMAND: $sff2fasta\n";
  system($sff2fasta);
  # Is not necessary because we have to remove the N's at the end of the reads from the fasta file.
  # my $sff2qual="sffinfo -q $outputTRIM > $qual";
  # print "COMMAND: $sff2qual\n";
  # system($sff2qual);
  # Remove N's at the end of the reads from fasta File
  $ind = rindex($fasta, '.');
  ## Create a temporal file with the sequence in one line
  my $onelineFasta=substr($fasta, 0, $ind).".oneline";
  print "Temp fasta file: $onelineFasta\n";
  open($fh, '>', $onelineFasta) or die "Could not open file '$onelineFasta' $!";
  open(IN,"<", $fasta) || die ("Error opening $fasta $!");
  my $firstLine=0;
  while (my $line = <IN>){
    chomp $line;
    if ($line=~/^>/) {
      if ($firstLine==0){
        print $fh $line."\n";
        $firstLine=1;
        }
      else{
        print $fh "\n".$line."\n";
        }
      }
    else{
      print $fh $line;
      }
    }
  print $fh "\n";
  close($fh);
  my $NtrimPoints=substr($fasta, 0, $ind).".NTP"; ## NTP: N-TrimPoints
  open($fh, '>', $NtrimPoints) or die "Could not open file '$NtrimPoints' $!";
  open (FH, "<", $onelineFasta) or die "Could not open file '$onelineFasta' $!";
  open (TH, "<", $finalTP) or die "Could not open file '$finalTP' $!";
  my $numN=0;
  my $linesSeq=0;
  my $linesTrim=0;
  while (my $header = <FH>){
    my $sequence = <FH>;
    chomp $sequence;
    $numN=0;
    $linesSeq++;
    my $revSeq=reverse $sequence;
    my @char = split //, $revSeq;
    foreach (@char){
      my $c=$_;
      if ($c ne "N"){
        last;
        }
      else{
        $numN++;
        }
      }
    my $ending=$#char-$numN+1;
    my $trim=<TH>;
    $linesTrim++;
    my @array = split(/\s+/, $trim);
    my $readID=$array[0];
    my $ini=$array[1];
    my $end=$ending+$ini-1;
    print $fh "$readID\t$ini\t$end\n";
    }
  close(FH);
  close(TH);
  close($fh);
  #--- Apply the new Trim Points
  my $lastTRIM=substr($fasta, 0, $ind)."-Ns";
  my $outputSFF=substr($fasta, 0, $ind);
  my $getLastReadsTRIM="sfffile -o $lastTRIM -tr $NtrimPoints $outputSFF";
  print "COMMAND: $getLastReadsTRIM \n";
  system($getLastReadsTRIM);
  my $lastFasta=substr($fasta, 0, $ind)."-Ns.fasta";
  my $getLastFasta="sffinfo -s $lastTRIM > $lastFasta";
  print "COMMAND: $getLastFasta \n";
  system($getLastFasta);
  my $outputQual=substr($fasta, 0, $ind)."-Ns.qual";
  my $getQual="sffinfo -q $lastTRIM > $outputQual";
  print "COMMAND: $getQual \n";
  system($getQual);
  }
print "DONE\n";
rmtree "cleaning_1";
foreach my $f (glob "err_seqcl_*"){
  unlink $f or print "unable to delete $f\n";
  }
#foreach my $f (glob "seqcl_*"){
 # unlink $f or print "unable to delete $f\n";
 # }
foreach my $f (glob "*.cidx"){
  unlink $f or print "unable to delete $f\n";
  }
foreach my $f (glob "*.sort"){
  unlink $f or print "unable to delete $f\n";
  }
#unlink $onelineFasta;
exit (0);
