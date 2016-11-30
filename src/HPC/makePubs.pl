#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs

my $seriesName = "CD3ZetaElectrostaticsRwall6"; # command line

my $i0Max = 13;
my $i1Max = 13;

for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
	my $WELLDEPTH = 10**(-4+0.5*($i0-1));

	for (my $i1 = 1; $i1 <= $i1Max; $i1++)
	{
		my $DEBYE = 10**(-4+0.5*($i1-1));

        my $fileName = $seriesName . "WellDepth" . "." . $i0 . "Debye" . "." . $i1;
		my $runName =  $fileName . ".pub";
	
		open (FOOD, ">pubs/$runName" );
		print FOOD << "EOF";

#!/bin/bash
#\$ -N $seriesName
#\$ -q free40i,free32i,free64
#\$ -e logs/$runName.err
#\$ -o logs/$runName.log

cd /pub/lclemens/polymer-c_runs/Dec012016ElectrostaticsRwall6Sweep/$seriesName

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Run your executable and then stuff
./metropolis.out parameters.txt $fileName -1 -1 $WELLDEPTH $DEBYE

echo Finished at `date`

EOF
		close FOOD;	
	}	
}
