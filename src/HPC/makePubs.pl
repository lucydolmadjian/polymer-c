#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs

my $seriesName = shift; # command line
my $N          = shift; 

my $i0Max = 23;
my $i1Max = 23;

for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
	my $c0 = 10**(-2+0.25*($i0-1));

	for (my $i1 = 1; $i1 <= $i1Max; $i1++)
	{
		my $c1 = 10**(-2+0.25*($i1-1));

		my $runName = $seriesName . "." . $N . "." . $i0 . "." . $i1 . ".pub";
	
		open (FOOD, ">pubs/$runName" );
		print FOOD << "EOF";

#!/bin/bash
#\$ -N $seriesName
#\$ -q free40i,free32i,free64
#\$ -e logs/$runName.err
#\$ -o logs/$runName.log

cd /pub/jallard/science/projects/LSM/1482/polymer/runs/$seriesName

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Run your executable and then stuff
./driveMetropolis params$N.txt $runName 0 $c0 $c1

echo Finished at `date`

EOF
		close FOOD;	
	}	
}
