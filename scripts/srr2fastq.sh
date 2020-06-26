#!/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/bash

for f in SRR1*; 
do 
	fasterq-dump --skip-technical --split-files $f; 
done
