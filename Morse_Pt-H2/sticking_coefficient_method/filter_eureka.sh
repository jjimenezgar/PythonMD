#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

# Nombre del archivo
file="best_parameters.txt"

# Saltar la primera línea (encabezado)
tail -n +2 "$file" | while IFS="," read EA Energy De1 Alfa De2 Alfa2 Re2; do
    DE="$De1"
    DE2="$De2"
    alfa="$Alfa"
    alfa2="$Alfa2"
    re2="$Re2"


cat ekalfa_*_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt > c1_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt 
cat probability_*_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt > c2_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt

paste c2_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt c1_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt > ek_sc_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt

# Asumimos que la columna con los números a comparar se encuentra en la primera posición
# del archivo

# Asumimos que la columna con los números a comparar se encuentra en la primera posición
# del archivo

# Obtenemos la primera línea del archivo
linea=$(head -n 1 ek_sc_"$DE"_"$DE2"_"$alfa"_"$alfa2"_"$re2".txt)

# Obtenemos el primer número de la línea
numero=$(echo $linea | awk '{print $1 + 0}')

if [ $(echo "$numero > 0.05" | bc) -eq 1 ]; then

# Si el número es menor que 0.15, continuamos con el siguiente ciclo
if [ $(echo "$numero < 0.15" | bc) -eq 1 ]; then

# Inicializamos la variable "mayor" con el primer número de la columna
mayor=$numero

# Iteramos sobre cada línea del archivo
while read linea; do
  # Obtenemos el primer número de la línea actual
  numero=$(echo $linea | awk '{print $1 + 0}')

  # Si el número de la línea actual es mayor que el número mayor actual,
  # actualizamos el valor de la variable "mayor"
if [ $(echo "$numero >= $mayor" | bc) -eq 1 ]; then
    mayor=$numero
	ek=$(echo $linea | awk '{print $2}')

fi

done < ek_sc_"$DE"_"$DE2"_"$alfa"_""$alfa2"_$re2".txt


if [ $(echo "$ek >= 60" | bc) -eq 1 ]; then
  echo "$mayor $ek $DE $DE2 $alfa $alfa2 $re2 " >> Eureka_1.txt 
fi
fi
fi

# Imprimimos el número mayor
done







