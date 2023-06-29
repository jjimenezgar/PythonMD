from calendar import c
from itertools import count
import numpy as np
import MDAnalysis as mda
import numpy.linalg
from statistics import mean
import sys


def atom_contact(file,atom1,atom2,r_cut):

    '''
    Esta función calcula el número de atomos que estan en contacto con la superficie"
    file: Archivo de entrada
    atom1: Atomo de referencia
    atom2: Atomo em contacto
    r_cut: radio de corte

    '''

    u=mda.Universe(file)

    term1 = u.select_atoms('name {name_atom1}'.format(name_atom1=atom1))
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2))

    number_atom1= len(term1)
    number_atom2=len(term2)


    list_contact=[]
    for ts in u.trajectory: ## lee cada frame
        count=0  # Cuenta el número de atomos en contacto
        for i in range(number_atom1): ## Empieza por comparar el primer atomo con el resto
            list_d=[]
            for x in range(number_atom2): ### Compara el atomo anterior con la seleccción 2
                r=term2[x].position-term1[i].position  # Calcula la diferencia vectorial
                d=numpy.linalg.norm(r) #Calcula la distancia
                if d < float(r_cut): # si la distancia es menor a X
                    list_d.append(d)  ## añade los atomos que cuenta los coloca en una lista
        
        

            if  len(list_d) > 0 and min(list_d) < float(r_cut):  ## toma la menor distancia de la lista de atomos que cuenta anteriormente 
                    count += 1  ## Cuenta los atomos en contacto
        
        list_contact.append(count) ## Crea una lista con los atomos en contacto

    return (mean(list_contact))

def cm_contact(file,atom1,atom2,r_cut):

    '''
    Esta función calcula el número de atomos que estan en contacto con la superficie"
    file: Archivo de entrada
    atom1: Atomo de referencia
    atom2: Atomo em contacto
    r_cut: radio de corte
    '''

    u=mda.Universe(file)

    term1 = u.select_atoms('index {name_atom1}'.format(name_atom1=atom1)) # Seleccionar grupo de atomos para calcular centro de masa
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2)) # Seleccionar superficie
    
    # MDAnalysis asigna las masas de los atomos automaticamente si estan bien escritos los elementos

    number_atom2=len(term2)


    list_contact=[]
    for ts in u.trajectory: ## lee cada frame
        count=0  # Cuenta el número de atomos en contacto
        list_d=[]
        cm = term1.center_of_mass() ## Centro de masa
        for x in range(number_atom2): ### Compara el centro de masa con la superficie
                r=cm-term2[x].position  # Calcula la diferencia vectorial
                d=numpy.linalg.norm(r) #Calcula la distancia        
                if d < float(r_cut): # si la distancia es menor a X
                    list_d.append(d)  ## añade los atomos que cuenta los coloca en una lista
        if  len(list_d) > 0 and min(list_d) < float(r_cut):  ## toma la menor distancia de la lista de atomos que cuenta anteriormente 
            count += 1  ## Cuenta los atomos en contacto
        list_contact.append(count) ## Crea una lista con los atomos en contacto

    splits = np.array_split(list_contact,50)  ### variable seria el numero de veces que se repite el calculo

    touch_count=0
    for array in splits:
        if sum(list(array)) > 10:
            touch_count += 1  ###  Si la molécula esta cerca mas de 5 frames(0.12 ps), se considera en contacto

    probability=round(touch_count/50, 4)  ## tomamos la probabilidada de cada trayectoria en contacto y la dividimos por el total

    return probability


def surface_contact(file,atom1,atom2,r_cut):

    '''
    Esta función calcula el número de atomos que estan en contacto en la superficie"
    file: Archivo de entrada
    atom1: Atomo de referencia
    atom2: Atomo em contacto
    r_cut: radio de corte
    '''

    u=mda.Universe(file)

    term1 = u.select_atoms('name {name_atom1}'.format(name_atom1=atom1)) # Seleccionar la superficie
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2)) # Seleccionar grupo de atomos para calcular centro de masa
    
    # MDAnalysis asigna las masas de los atomos automaticamente si estan bien escritos los elementos

    number_atom1= len(term1)
    number_atom2=len(term2)

   
    list_contact=[]
    for ts in u.trajectory:
        count=0
        for x in range(number_atom1):  ### Cambia el atomo a evaluar
            list_d=[]
            for i in range(number_atom2):
                r=term1[x].position-term2[i].position
                d=numpy.linalg.norm(r)
                if d < float(r_cut):
                    list_d.append(d)
        
        

            if len(list_d) > 0 and   min(list_d) < float(r_cut):
                    count += 1
        
        list_contact.append(count)

    return (mean(list_contact))



if __name__ == '__main__':

    '''
    Esta función calcula el número de atomos que estan en contacto en la superficie"
    file: Archivo de entrada
    atom1: Atomo de referencia
    atom2: Atomo em contacto
    r_cut: radio de corte
    '''

    r=(globals()[sys.argv[1]](sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]))   ### abrir función desde la terminal
    print(r)
