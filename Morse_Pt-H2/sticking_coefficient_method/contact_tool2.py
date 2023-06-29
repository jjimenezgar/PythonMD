from calendar import c
from itertools import count
import numpy as np
import MDAnalysis as mda
import numpy.linalg
from statistics import mean
import sys


def atom_contact(file, atom1, atom2, r_cut):
    '''
    This function calculates the number of atoms in contact with the surface.
    
    Parameters:
    file: Input file
    atom1: Reference atom
    atom2: Atom in contact
    r_cut: Cutoff radius
    
    Returns:
    The mean number of atoms in contact with the surface.
    '''
    u = mda.Universe(file)

    term1 = u.select_atoms('name {name_atom1}'.format(name_atom1=atom1))
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2))

    number_atom1 = len(term1)
    number_atom2 = len(term2)

    list_contact = []
    for ts in u.trajectory:
        count = 0
        for i in range(number_atom1):
            list_d = []
            for x in range(number_atom2):
                r = term2[x].position - term1[i].position
                d = numpy.linalg.norm(r)
                if d < float(r_cut):
                    list_d.append(d)
        
            if len(list_d) > 0 and min(list_d) < float(r_cut):
                count += 1
        
        list_contact.append(count)

    return mean(list_contact)


def cm_contact(file, atom1, atom2, r_cut):
    '''
    This function calculates the number of atoms in contact with the surface.
    
    Parameters:
    file: Input file
    atom1: Reference atom
    atom2: Atom in contact
    r_cut: Cutoff radius
    
    Returns:
    The contact probability based on the center of mass of atom1 with the surface.
    '''
    u = mda.Universe(file)

    term1 = u.select_atoms('index {name_atom1}'.format(name_atom1=atom1))
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2))

    number_atom2 = len(term2)

    list_contact = []
    for ts in u.trajectory:
        count = 0
        list_d = []
        cm = term1.center_of_mass()
        for x in range(number_atom2):
            r = cm - term2[x].position
            d = numpy.linalg.norm(r)
            if d < float(r_cut):
                list_d.append(d)
        
        if len(list_d) > 0 and min(list_d) < float(r_cut):
            count += 1
        
        list_contact.append(count)

    splits = np.array_split(list_contact, 50)
    touch_count = 0
    for array in splits:
        if sum(list(array)) > 10:
            touch_count += 1

    probability = round(touch_count / 50, 4)

    return probability


def surface_contact(file, atom1, atom2, r_cut):
    '''
    This function calculates the number of atoms in contact with the surface.
    
    Parameters:
    file: Input file
    atom1: Reference atom
    atom2: Atom in contact
    r_cut: Cutoff radius
    
    Returns:
    The mean number of atoms in contact on the surface.
    '''
    u = mda.Universe(file)

    term1 = u.select_atoms('name {name_atom1}'.format(name_atom1=atom1))
    term2 = u.select_atoms('name {name_atom2}'.format(name_atom2=atom2))

    number_atom1 = len(term1)
    number_atom2 = len(term2)

    list_contact = []
    for ts in u.trajectory:
        count = 0
        for x in range(number_atom1):
            list_d = []
            for i in range(number_atom2):
                r = term1[x].position - term2[i].position
                d = numpy.linalg.norm(r)
                if d < float(r_cut):
                    list_d.append(d)
        
            if len(list_d) > 0 and min(list_d) < float(r_cut):
                count += 1
        
        list_contact.append(count)

    return mean(list_contact)


if __name__ == '__main__':
    '''
    This function calculates the number of atoms in contact with the surface.
    
    Parameters:
    file: Input file
    atom1: Reference atom
    atom2: Atom in contact
    r_cut: Cutoff radius
    '''
    r = (globals()[sys.argv[1]](sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]))
    print(r)
