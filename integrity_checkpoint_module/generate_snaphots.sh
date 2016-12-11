#!/bin/bash

generate_snapsot()
{
    rank=$1
    phase=0

    mkdir $rank

    cd $rank

    for index in `seq 0 9`;
    do
        echo $phase"_"$index
        touch $phase"_"$index
    done

    cd ..
}

########################## main ##############################################

rank=$1

echo "Snapshot name: [checkpoint_phase]_[checkpoint_index]"
echo "Generate snaphots for ranks form 0 to $rank"

mkdir "snapshot"

cd "snapshot/"

for i in `seq 0 $rank`;
do
    echo ""
    echo "Generate snaphots for rank $i"
    generate_snapsot $i
done

