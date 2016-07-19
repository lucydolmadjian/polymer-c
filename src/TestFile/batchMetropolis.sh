# Parallelization Script
# August 11, 2015

NRequested=0            # initialize number of runs submitted

NRODS=143

RATIO=10

STIFFENRANGE=3

ITERATIONS=1

TOTALITERATIONS=`wc -l < PhosphorylatediSitesNoSpace.txt`

echo "Length of file is $TOTALITERATIONS"



    # set number of runs submitted by checking the running processes for lines with the program name, set NRequested to number of lines (may include one more than actual number of runs, since it counts grep -c metropolis in its tally)

#NRequested=`ps | grep -c metropolis`

# while number of iterations ran is less than or equal to total number of iterations desired, loop through runs

while (( $ITERATIONS <= $TOTALITERATIONS ))
    do

    # loop to periodically check how many runs are submitted

#while (( $NRequested >= $1 ))   # while number requested is greater than number of processors we want to use (user input in command line)

#do
#sleep 1     # wait 1 sec before checking again

#NRequested=`ps | grep -c metropolis` # check again

#done

    # loop to submit more runs until reach max number of processors we want to use ($1)

#while (( $NRequested < $1 && $ITERATIONS <= $TOTALITERATIONS ))
#do

            # read specified line of text file

            STIFFISITES=`awk 'NR==iter' iter=$ITERATIONS PhosphorylatediSites.txt`

            # print to screen the line read
            echo "Line $ITERATIONS of file is $STIFFISITES"

#./raninit.out       # initialize random number generator

            # run program with specified parameters
#./metropolis.out StiffenedSegments.$ITERATIONS $NRODS $RATIO $STIFFISITES $STIFFENRANGE &

            # If user gives V or v as second command line argument, then code will be verbose. Any other input will result in non-verbose.
#if [[ $2 == "V" || $2 == "v" ]]
#then
                    # print to screen the process ID and the name of the run
#echo "PID of StiffenedSegments.$ITERATIONS is $!"
#fi

            # update number of running programs
#NRequested=`ps | grep -c metropolis`

            # increase iteration run by 1
ITERATIONS=$(( ITERATIONS + 1 ))
#done

done

# wait for all background processes to finish before concatenating files
#wait

# loop through all files, concatenate them into one file
# for ((N=1; N<=$NRODSMAX; N++))
# do
# for ((R=1; R<=$RATIOMAX; R++))
# do

#cat StiffenedSegments.* >> StiffenedSegmentsSweep.txt

# done
# done