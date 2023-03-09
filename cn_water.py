#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:17:47 2023

@author: johnconnor
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def cn_atom(file,r_cut_lo,r_cut_hi):
# Load trajectory and topology
    u = mda.Universe(file)

# Select only oxygen atoms
    oxygen = u.select_atoms("prop 5 < z and prop z < 30 and name O1 and around {} name Cgr and not around {} name Cgr".format(r_cut_hi,r_cut_lo))
    
# Loop over all frames and calculate tetrahedral order parameter

    coord=[]
    for ts in u.trajectory:
        if ts.frame == 0:
            for atom in oxygen:
                neighbors = u.select_atoms("name O1 and around 3.5 index {}".format(atom.index))        
                coord.append(len(neighbors))
    
    return np.mean(coord)

def graf_cn(X,Y):
    
    plt.plot(X, Y, color='red', linestyle='dashed', linewidth=2)

   
    plt.title('Water molecules CN in the pore', fontsize=20)

    # Agregar etiquetas a los ejes
    plt.xlabel('Radius (Å)')
    plt.ylabel('Cordination Number (CN)')

    # Mostrar el grafico
    plt.show()
    
    

def processing(olig):
    file = "/home/johnconnor/Documentos/Mesoporos/cluster_result/OMC3/{}_300.xyz".format(olig)
#    image_name="{}_omc3_pore_msd_300.png".format(olig)
    df = pd.DataFrame(columns=["R", "CN"])

   
    radio_p=15
    step=3
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
    
       for olig in [4,8,12,24,36]:
          #olig="water"
          processing(olig)
        
                   
        


        
    