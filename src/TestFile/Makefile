all:
	gcc -O3 driveMetropolis.c -o metropolis.out -lm

	chmod a+x batchMetropolis_Multiple.sh

	cp metropolis.out ~/Documents/polymer-c_runs/.

multiple:

	gcc -O3 driveMetropolis.c -o metropolis.x -lm

	chmod a+x batchMetropolis_Multiple.sh

	./batchMetropolis_Multiple.sh 7 v

stiffen:

	gcc -O3 driveMetropolis.c -o metropolis.out -lm

	chmod a+x batchMetropolis_Stiffen.sh

	./batchMetropolis_Stiffen.sh 7 v

