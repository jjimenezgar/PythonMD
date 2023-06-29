for i in {7..10}; do

./stick_new.sh

cat probability_* > probability.txt

cat ekalfa_* > ek_alfa.txt

paste ek_alfa.txt probability.txt > sco_$i.txt

rm probability*
rm ekalfa_*
rm ek_alfa.txt

done

cat sco_* > scofinal.txt

rm sco_*
