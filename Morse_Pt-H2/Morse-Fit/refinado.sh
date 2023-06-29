#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

for DE in $(seq 0.201 0.01 0.40); do
    for DE2 in 8e-7 9e-7 1e-6 1.1e-6 1.2e-6; do
	for alpha in $(seq 4.01 0.2 5.0); do
	   for alpha2 in $(seq 0.301 0.02 1.0); do
		for re2 in $(seq 10.01 0.2 20); do		

awk '$1 >= 2.0 && $1 <= 3.5 {print $1, $2}' Energy_"$DE"__"$alpha"__2.2__"$DE2"__"$alpha2"__"$re2"_.dat > list_ea_"$DE"_"$alpha"_"$DE2"_"$alpha2"_"$re2".dat
awk '$1 >= 1.5 && $1 <= 2.0 {print $1, $2}' Energy_"$DE"__"$alpha"__2.2__"$DE2"__"$alpha2"__"$re2"_.dat > list_e_"$DE"_"$alpha"_"$DE2"_"$alpha2"_"$re2".dat


Ea=$(awk 'BEGIN{a=0}{if ($2>0+a) a=$2} END{print a}' list_ea_"$DE"_"$alpha"_"$DE2"_"$alpha2"_"$re2".dat)	
Energy=$(awk 'BEGIN{a=1000}{if ($2<0+a) a=$2} END{print a}' list_e_"$DE"_"$alpha"_"$DE2"_"$alpha2"_"$re2".dat)

echo "$Ea $Energy $DE $alpha $DE2 $alpha2 $re2" >> parameters_morse.txt

done
done
done
done
done
