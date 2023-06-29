#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

# File name
file="best_parameters.txt"

# Skip the first line (header)
tail -n +2 "$file" | while IFS="," read EA Energy De1 Alfa De2 Alfa2 Re2; do
    DE="$De1"
    DE2="$De2"
    alpha="$Alfa"
    alpha2="$Alfa2"
    re2="$Re2"

    for ek in 33 44 54 63; do
        sed -e 's/EK/'$ek'/g' in.Pt_H | sed -e 's/DE/'$DE'/g' | sed -e 's/D2/'$DE2'/g' | sed -e 's/alfa/'$alpha'/g'| sed -e 's/alpha2/'$alpha2'/g'  | sed -e 's/RE/'$re2'/g'  > in.Pt_H_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2"

        for test in {1..50}; do
            R=$(($RANDOM%180+1))
            xran=$(($RANDOM%20+2))
            yran=$(($RANDOM%20+2))

            rad=$(bc -l <<< "scale=4; $R*(3.14/180)")
            x=$(bc -l <<< "scale=4;$xran + 0.75* c($rad)*c($rad)")
            y=$(bc -l <<< "scale=4;$yran + 0.75* c($rad)*s($rad)")
            z=$(bc -l <<< "scale=4;23.30 + 0.75*s($rad)")

            sed s/coord_x/$x/g data_Pt_H | sed s/x_ran/$xran/g | sed s/y_ran/$yran/g | sed s/coord_y/$y/g | sed s/coord_z/$z/g > update.data

            lmp_mpi <  in.Pt_H_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2" -v rnseed $test > in.Pt_H_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2".out
	
	        echo "$ek" "$DE" "$DE2" "$alpha" "$alpha2" "$re2". > ekalfa_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2".txt
        done

        python /home/contact_tool2.py cm_contact SC_Pt_H_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2".xyz 601:602 Pt 3 > probability_"$ek"_"$DE"_"$DE2"_"$alpha"_"$alpha2"_"$re2".txt 

        rm SC_Pt_H_* in.Pt_H_*
    done
done
