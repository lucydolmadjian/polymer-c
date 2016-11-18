#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs

my $seriesName = shift; # command line
my $N          = shift; 

my $i0Max = 23;
my $i1Max = 23;

for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
	my $i1 = 1;
	while ( $i1 <= $i1Max ) 
	{

			my $runName = $seriesName . "." . $N . "." . $i0 . "." . $i1 . ".pub";
			system "qsub pubs/$runName";

			print $i1;
			$i1++;
			sleep(3);	

	}	
}
