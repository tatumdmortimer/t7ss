#!/bin/bash

/opt/PepPrograms/prokka-1.11/bin/prokka --prefix $1_$2_$3 --locustag $1_$2_$3 --genus $1 --species $2 --strain $3 --usegenus $1_$2_$3.fasta
tar -czvf $1_$2_$3.tar.gz $1_$2_$3/
