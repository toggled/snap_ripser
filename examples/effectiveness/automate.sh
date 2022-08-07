#!/bin/sh

#  automate_timing.sh
#  
#
#  Created by Naheed Anjum Arafat on 6/10/19.
#
BASEDIR=$(dirname $0)
datadir=$BASEDIR"/../"
declare -a datanames
# datanames=("celegans_final.edges" "power.edge")
datanames=("celegans_final.edges")
declare -a diams
diams=(1.33 46)
declare -a V
V=(4941 317080 334863)
declare -a algs
algs=("epsnetring" "eps_filtering" "eps_baseline2")
declare -a mn
#mn=(4 5 10)
mn=(0.01 1)
inc=(0.05 1)

### Run eps-ring first
#echo "Running eps-ring\n"
i=0
for f_name in "${datanames[@]}"; do
    echo "$datadir$f_name";
    for k in {0..20}; do
        lower_lt=$(bc <<< "scale=3; ${mn[$i]}+$k*${inc[$i]}")
        echo $lower_lt
        for algo in "${algs[@]}"; do
            if [ "$algo" == "${algs[0]}" ]; then
#               echo $algo
                ./lazywitnessph_effect --format wgraph "$datadir$f_name" --heuristic $algo --dim 1 --epsilon $lower_lt --iterations 20
    #            done
            fi
            if [ "$algo" == "${algs[1]}" ]; then
#                echo $algo
#                echo $lower_lt
                ./lazywitnessph_effect --format graph "$datadir$f_name" --heuristic $algo --dim 1 --epsilon $lower_lt --iterations 20
            fi
        done
    done
    i=$((i+1))
done

#echo "Running eps_filtering\n"
#i=0
#for f_name in "${datanames[@]}"; do
#echo "$datadir$f_name";
#for k in {0..20}; do
#lower_lt=$((diams[i]/10 + k))
##        echo "$((lower_lt)) "$((diams[i]/10))
#./lazywitnessph_effic --format graph "$datadir$f_name" --heuristic eps_filtering --dim 1 --epsilon $lower_lt
#done
#i=$((i+1))
#done
#
#echo "Running eps_kruskal\n"
#i=0
#for f_name in "${datanames[@]}"; do
#echo "$datadir$f_name";
#for k in {0..20}; do
#lower_lt=$((diams[i]/10 + k))
##        echo "$((lower_lt)) "$((diams[i]/10))
#./lazywitnessph_effic --format graph "$datadir$f_name" --heuristic eps_kruskal --dim 1 --epsilon $lower_lt
#done
#i=$((i+1))
#done
#
#echo "Running eps_baseline\n"
#i=0
#for f_name in "${datanames[@]}"; do
#echo "$datadir$f_name";
#for k in {0..20}; do
#lower_lt=$((diams[i]/10 + k))
##        echo "$((lower_lt)) "$((diams[i]/10))
#./lazywitnessph_effic --format graph "$datadir$f_name" --heuristic epsbaseline2 --dim 1 --epsilon $lower_lt
#done
#i=$((i+1))
#done
