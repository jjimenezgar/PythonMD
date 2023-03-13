#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Este script calcula el número de coordinación de un conjunto de átomos de oxígeno
en una simulación de dinámica molecular utilizando la librería MDAnalysis.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def cn_atom(file, r_cut_lo, r_cut_hi):
    """
    Calcula el número de coordinación de un conjunto de átomos de oxígeno.

    Parameters
    ----------
    file : str
        Ruta del archivo que contiene la simulación.
    r_cut_lo : float
        Radio de corte inferior para buscar átomos de Pt.
    r_cut_hi : float
        Radio de corte superior para buscar átomos de Pt.

    Returns
    -------
    float
        El número de coordinación medio del conjunto de átomos de oxígeno.

    """
    # Carga la trayectoria y la topología
    u = mda.Universe(file)

    # Selecciona solo los átomos de oxígeno
    oxygen = u.select_atoms("prop 30 < z and prop z < 70 and name O1 and around {} name Pt and not around {} name Pt".format(r_cut_hi, r_cut_lo))
    
    coord = [len(u.select_atoms("name O1 and around 3.25 index {}".format(atom.index))) 
         for ts in u.trajectory 
         if ts.frame == 0 
         for atom in oxygen]
         
    return np.mean(coord)



def graf_cn(X, Y):
    """
    Grafica el número de coordinación en función del radio.

    Parameters
    ----------
    X : array-like
        Arreglo con los radios.
    Y : array-like
        Arreglo con los números de coordinación correspondientes.

    Returns
    -------
    None.

    """
    plt.plot(X, Y, color='black', linestyle='dashed', linewidth=2, label="12")

    plt.title('Coordination number/water', fontsize=18)

    # Agrega etiquetas a los ejes
    plt.xlabel('Radius (Å)', fontsize=16)
    plt.ylabel('Coordination Number (CN)', fontsize=16)
    # plt.xticks(np.arange

    plt.legend(fontsize=12)

    # Mostrar el grafico
    plt.show()
    
    

def processing(olig):
    file = "/home/johnconnor/Documentos/Mesoporos/cluster_result/OMC3/{}_cn.xyz".format(olig)
    df = pd.DataFrame(columns=["R", "CN"])

   
    radio_p=20
    step=5
    final= radio_p + step
    
    for i in np.arange(step, final, step):
        # Calcular el valor de r_cut
        r_cut_lo = (final) - i
        r_cut_hi = radio_p - i


        
        CN = cn_atom(file,r_cut_hi,r_cut_lo)

    
        # Agregar una fila al DataFrame con los valores de r_cut y "difusión" actuales
        df = df.append({"R": r_cut_lo, "CN": CN}, ignore_index=True)
        
    graf_cn(df.R,df.CN)
    
if __name__ == '__main__':
    
       for olig in [12]:
          #olig="water"
          processing(olig)
        
                   
        


        
    