#!/bin/bash

set -e 

generate_sequence() {
    local len=$1
    local seq=""
    local bases=(A C G T)
    for (( i=0; i<$len; i++ )); do
        seq+=${bases[$RANDOM % 4]}
    done
    echo "$seq"
}


sequence0="GCAAGTGACAA"
sequence1="GCAAGTGACAATTCCTGAGAATAAT"
sequence2="TACGACAGCGAATTCCTGAGAATAATTCGAAAATTCGGACCTCGTTCGTCAATAGCTGTCATAAGGTGGGGTTACATCCCGCTCAGTCTATAATAGTGCACTTTTGTGCAGAGATTGCTGAGTGGCGAGAATTACTGTCTGGGAGCTAATGCTCACCGGGTTCAAGTGGACACAACTCGGATCCTAATACGCACATATACCGTAAATTCACGTCAACATCTCGGTATGCGTAGACACAATTGATCAATATCCTGAAAGCCCCG"
sequenceN="TACGACAGCGAATTCCTGAGAATAATTCGAAAATTCGGACCTCGTTCGTCAATAGCTGTCATAAGGTGGGGTTACATCCCGCTCAGTCTATAATAGTGCACTTTTGTGCAGAGATTGCTGAGTGGCGAGAATTACTGTCTGGGAGCTAATGCTCACCGGGTTCAAGTGGACACAACTCGGATCCTAATACGCACATATACCGTAAATTCACGTCAACATCTCGGTATGCGTAGACACAATTGATCAATATCCTGAAAGCCCCGTGATTAAGCTGCGGGGGAATGGCTTATACTCATTAGTGAATAAATACATAGGCGCCAGAGATTATGAACGTTCCTAAGAGTTGGATACCACCATTCAAGGGTTCACGCGCCGGTGTATTCGACTCATCTACGGCCATCAGTGGCGAGTTTACTTACGTGTTAACAGAGTACCGCCCGATTTTCCATGGGGAGTGTATTCAGATGATGCGGGAGACCGGGCAGTAAAATCGCCCCCATCTGAGAATGGCGATCCTTGTGCGTGTCGGTTCGCATTTTTGCTGAGACTAAAAAGACTCCAAAATGTAGATATTATGGTAGCTTGGGTGATGGCCAGTTTACATAGACAGTAAATTAAGTTACTCGAAACCGTACTCATTGTTGTGGAGCCGAC"

a=5
b=2


for i in {1..200} #10
do

    a=$(( RANDOM % 989 + 11 ))  # Random integer between 11 and 999 inclusive

    # Calculate 'b' constraints: 6 < b < a and b < 100
    b_min=2
    b_max=$(( a - 1 ))
    if (( b_max > 31 )); then
        b_max=31
    fi

    # Ensure b_min does not exceed b_max
    if (( b_max < b_min )); then
        b=$b_min
    else
        b=$(( RANDOM % (b_max - b_min + 1) + b_min ))
    fi
    sequence=$(generate_sequence 10000)
    # Execute the test command with the generated parameters
    echo "./test \"$sequence\" $a $b 0"
    ./test "$sequence" $a $b 0
    # ./bin/syncmer_tree ./bin/csyncmer "$seqfile" $a $b
done


# sequence=CACACAGA
# a=5
# b=2

# echo "./test \"$sequenceN\" $a $b"
# ./test "$sequenceN" $a $b 0

echo "All tests for correctness are OK. Now testing speed."

echo "Setting cpu to 2.6 GHz."
echo "sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz"
sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz

echo "Disabling HYPERTREADING for more accurate time estimation"
echo "sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'"
sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'

echo "RUNNING SPEED TEST"
echo "./test data/chr19_bit.fa 31 11 1"
./test data/chr19_bit.fa 31 11 1

echo "Re-enabling HYPERTREADING."
echo "sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'"
sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'