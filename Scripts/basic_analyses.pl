#!/usr/bin/env perl

############     basic_analyses.pl     ############
#
#	This Perl script performs some basic analyses:
#	- count number of sequences
#	- check sequence length distribution
#	- check sequence base composition
#	
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de



############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed and
#	the output file name as arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	The output file name has to be passed to the script via -o:
#	-o output_file.txt
#
#	If no output file name is passed to the script, the results
#	will be displayed in the command promt.
#
#	For example you can type the following command:
#	perl basic_analyses.pl -i file1.fas -o output_file.txt
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	Multiple files and combinations of all kinds of arguments are
#	allowed:
#	perl basic_analyses.pl -i file1.fas -I list_of_files.txt -o output_file.txt



@input_files=();
$|=1;

###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl basic_analysis.pl -argument1 -argument2 ...\n\n";
	exit;
	}

$argv="";
foreach(@ARGV)
	{
	$argv.=$_;
	}
@arguments=split('-',$argv);

foreach(@arguments)
	{
	if($_=~/^ *i/)
		{
		$_=~s/^ *i//;
		$_=~s/ //g;
		push(@input_files,$_);
		}
	elsif($_=~/^ *I/)
		{
		$_=~s/^ *I//;
		$_=~s/ //g;
		open(FILES_IN,"$_");
		while(<FILES_IN>)
			{
			unless($_=~/^\s*$/)
				{
				$_=~s/\s//sg;
				push(@input_files,$_);
				}
			}
		}
	elsif($_=~/^ *o/)
		{
		$_=~s/^ *o//;
		$_=~s/ //g;
		$output_file_name=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
	exit;
	}

###   PRINT ARGUMENTS   ###
print"The following files will be processed:\n";
foreach(@input_files)
	{
	if(-e $_)
		{
		print"$_\n";
		push(@input_files_ok,$_);
		}
	else
		{
		print"could not find file: $_. It will be ignored.\n";
		}
	}

###   START   ###
if($output_file_name)
	{
	open(OUT,">$output_file_name");
	}

foreach$file(@input_files_ok)
	{
	$number_of_sequences=0;
	$count_A=0;
	$count_T=0;
	$count_G=0;
	$count_C=0;
	$count_U=0;
	$count_N=0;
	$count_other=0;
	%length_distribution=0;
	
	print"processing $file\n\n";
	open(IN,$file);
	while(<IN>)
		{		
		if($_=~/^>/)
			{
			$number_of_sequences++;
			}
		elsif($_!~/^\s*$/)
			{
			chomp$_;
			$length_distribution{length$_}++;
			$count_A+=$_=~s/A//g;
			$count_T+=$_=~s/T//g;
			$count_G+=$_=~s/G//g;
			$count_C+=$_=~s/C//g;
			$count_U+=$_=~s/U//g;
			$count_N+=$_=~s/N//g;
			$count_other+=$_=~s/[\n]//g;
			}
		}
	close IN;
	
	$ATGCUN=$count_A+$count_T+$count_G+$count_C+$count_U+$count_N;
	$percent_A=$count_A/$ATGCUN*100;
	$percent_T=$count_T/$ATGCUN*100;
	$percent_G=$count_G/$ATGCUN*100;
	$percent_C=$count_C/$ATGCUN*100;
	$percent_U=$count_U/$ATGCUN*100;
	$percent_N=$count_N/$ATGCUN*100;

	@length_distribution=();
	foreach(keys(%length_distribution))
		{
		push(@length_distribution,$_);
		}
	@length_distribution=sort{$a<=>$b}@length_distribution;

	if($output_file_name)
		{
		print OUT "RESULTS $file\nNumber of sequences: $number_of_sequences\n\n";
		print OUT "Base\tabsolute\tpercent:\nA\t$count_A\t\t$percent_A\nT\t$count_T\t\t$percent_T\nG\t$count_G\t\t$percent_G\nC\t$count_C\t\t$percent_C\nU\t$count_U\t\t$percent_U\nN\t$count_N\t\t$percent_N\nother\t$count_other\n*other characters than ATGCUN are ignored for calculation of percentage.\n\n";
		print OUT "Length distribution (sequence length / number of sequences)\n";
		foreach(@length_distribution)
			{
			unless($_==0)
				{
				print OUT "$_\t$length_distribution{$_}\n";
				}
			}
		print OUT "\n-------------------------------------------\n";
		}
	else
		{
		print"RESULTS $file\nNumber of sequences: $number_of_sequences\n\n";
		print"Base\tabsolute\tpercent:\nA\t$count_A\t\t$percent_A\nT\t$count_T\t\t$percent_T\nG\t$count_G\t\t$percent_G\nC\t$count_C\t\t$percent_C\nU\t$count_U\t\t$percent_U\nN\t$count_N\t\t$percent_N\nother\t$count_other\n*other characters than ATGCUN are ignored for calculation of percentage.\n\n";
		print"Length distribution (sequence length / number of sequences)\n";
		foreach(@length_distribution)
			{
			unless($_==0)
				{
				print"$_\t$length_distribution{$_}\n";
				}
			}
		print"\n-------------------------------------------\n";
		}
	}
close OUT;
if($output_file_name)
	{
	print"Results saved at $output_file_name\n\n";
	}
exit;