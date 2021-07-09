#!/usr/bin/env python
'''
beta2omega.py: Backend to perform transformations from beta to omega

__author__     = 'Camilo A. F. Salvador'
__email__      = "cafsss@gmail.com"

# -- INFO -- #
The script works 3x3x3 bcc superell and direct coordinates (0-1);
It will paste the positions from a VASP input file ignoring l1-l9 e.g.:

l1    #Ti40Nb12Al2
l2    #3.33978
l3    #  -3.00000000  0.00000000  0.00000000
l4    #   0.00000000  0.00000000  3.00000000
l5    #   0.00000000  3.00000000  0.00000000
l6    #Ti  Nb  Al
l7    #40  12  2
l8    #Selective dynamics
l9    #Direct
ll    #  0.000000  0.166666  0.500000  T T T
ll    #  0.333333  0.666666  1.000000  T T T
ll    #  (...)

The header of the POSCAR will be pasted into the OUTCAR file;
The user selects the desired variant (1-4) with the folowing functions:
   - beta_to_omega1
   - beta_to_omega2
   - beta_to_omega3
   - beta_to_omega4

# -- USAGE -- #
import beta2omega
from beta2omega import *

print_to_file(beta_to_omega1(load_atoms()))
'''

# -- start --

import os
import sys
import numpy as np

script_dir  = os.path.dirname(__file__)
input_path  = os.path.join(script_dir, "POSCAR")
output_path = os.path.join(script_dir, "OUTCAR")

beta_pos = np.array([0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]) # the beta original positions
factor = 1  # full collapse; factor = 2, half collapse; 
delta = ((((1/6)/3)/2)/factor) # the expected +/- shift from original positions (per axis)

header = [""]  # header of strings
atoms  = []    # a list of objs (Atom)

class Atom(object):
	# Atom class allows storing the positions (x, y, z) of each atom;
    def __init__(self):
        self.number = int(420)
        self.pos = [0.0, 0.0, 0.0]
        self.deg = 0.0
        self.shift = 0.0

    def get_deg(self):
        return self.deg

    def __repr__(self):
    	# exact form to be printed in the OUTCAR file
    	vasp_output = "  {:.6f}  {:.6f}  {:.6f}  T T T\n".format(self.pos[0], self.pos[1], self.pos[2])
    	return vasp_output

def load_atoms():
    i = 1
    for line in open(input_path):
        if i < 10:
            header.append(line)
        else:
            atom = Atom()
            atom.number = i-9
            atom.pos[0] = float(line.split('  ')[1])
            atom.pos[1] = float(line.split('  ')[2])
            atom.pos[2] = float(line.split('  ')[3])
            atoms.append(atom)
        i += 1
    return header, atoms
		
def scale(a):
    return (1/6)*a

def shiftm1_abc(atom):
    # determine the shift multiplier (m) based on a+b+c
    # see @ref1
    #Invar       02,05,08,11,14
    #Up  -->  00,03,06,09,12,15,18
    #Down <-     04,07,10,13,16
    # this function generates the trivial omega variant #1
    m = 0
    a_neg_b = 6*(atom.pos[0] - atom.pos[1]) #rescale to compare 6*

    for c in [2,5,8,11,14]:
        if np.isclose(c, (sum(atom.pos))*6):
            m =  0
    
    for c in [0,3,6,9,12,15,18]:
        if np.isclose(c, (sum(atom.pos))*6):
            m = +1
    
    for c in [4,7,10,13,16]:
        if np.isclose(c, (sum(atom.pos))*6):
            m = -1
    
    return m 

def beta_to_omega1(atoms):
    # receives a list of atoms
    # returns the list with positions shifted
    for atom in atoms:
        atom.shift = delta*shiftm1_abc(atom)
        atom.pos[0] += atom.shift
        atom.pos[1] += atom.shift
        atom.pos[2] += atom.shift
    return atoms

def shiftm2_abc(atom):
    # determine the shift multiplier (m) based on a+b+c
    # see @ref1
    #Invar       04,07,10,13,16
    #Up  -->     02,05,08,11,14
    #Down <-  00,03,06,09,12,15,18
    # this function generates the trivial omega variant #1
    m = 0
    a_neg_b = 6*(atom.pos[0] - atom.pos[1]) #rescale to compare 6*

    for c in [4,7,10,13,16]:
        if np.isclose(c, (sum(atom.pos))*6):
            m =  0
    
    for c in [2,5,8,11,14]:
        if np.isclose(c, (sum(atom.pos))*6):
            m = +1
    
    for c in [0,3,6,9,12,15,18]:
        if np.isclose(c, (sum(atom.pos))*6):
            m = -1
    
    return m 

def beta_to_omega2(atoms):
    # receives a list of atoms
    # returns the list with positions shifted
    for atom in atoms:
        atom.shift = delta*shiftm2_abc(atom)
        atom.pos[0] += atom.shift
        atom.pos[1] += atom.shift
        atom.pos[2] += atom.shift
    return atoms

def shiftm3_abc(atom):
    # determine the shift multiplier (m) based on atom coordinates (abc)
    # see @ref2
    # Invar       222,555
    # Up  --> 000,333,666
    # Down <-     111,444
    # this function generates the omega variant #3
    m = 0
    a_neg_b = 6*(atom.pos[0] - atom.pos[1]) #rescale to compare 6*

    # If c=(2, 5), atoms are naturally invariant
    for c in [2, 5]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m =  0
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m = +1
            else: 
                m = -1

    # If c=(0, 3, 6), atoms are +shifted
    for c in [0, 3, 6]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m = +1
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m = -1
            else: 
                m =  0

    # If c=(1, 4), atoms are -shifted
    for c in [1, 4]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m = -1
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m =  0
            else: 
                m = +1
    return m

def beta_to_omega3(atoms):
    # receives a list of atoms
    # returns the list with positions shifted
    for atom in atoms:
        atom.shift = delta*shiftm3_abc(atom)
        atom.pos[0] += atom.shift
        atom.pos[1] -= atom.shift # more details below @ref3
        atom.pos[2] += atom.shift
    return atoms

def shiftm4_abc(atom):
    # determine the shift multiplier (m) based on atom coordinates (abc)
    # see @ref2
    # Invar       111,444
    # Down <- 000,333,666
    # Up  -->     222,555
    # this function generates the omega variant #4
    m = 0
    a_neg_b = 6*(atom.pos[0] - atom.pos[1]) #rescale to compare 6*

    # If c=(1, 4), atoms are naturally invariant
    for c in [1, 4]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m =  0
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m = +1
            else: 
                m = -1

    # If c=(0, 3, 6), atoms are -shifted
    for c in [0, 3, 6]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m = -1
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m =  0
            else: 
                m = +1

    # If c=(2, 5), atoms are +shifted
    for c in [2, 5]:
        if np.isclose(atom.pos[2], scale(c)):
            if np.isclose(a_neg_b, 0):
                m = +1
            # unless ...
            elif np.isclose(a_neg_b, -2) or np.isclose(a_neg_b, +4):
                m = -1
            else: 
                m =  0
    return m

def beta_to_omega4(atoms):
    # receives a list of atoms
    # returns the list with positions shifted
    for atom in atoms:
        atom.shift = delta*shiftm4_abc(atom)
        atom.pos[0] += atom.shift
        atom.pos[1] -= atom.shift # more details below @ref3
        atom.pos[2] += atom.shift
    return atoms


def print_to_file(header, atoms):
	# helper function to print atoms to a file
	outfile = open(output_path, "w")
	outfile.write(''.join(header))
	for atom in atoms:
		outfile.write(repr(atom))
	outfile.close()

# -- end --
