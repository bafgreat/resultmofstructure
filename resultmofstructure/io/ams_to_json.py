#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import glob
from ase.io import read
from ase import Atoms
import ams_to_ase
import numpy as np
import filetyper
import coordlibrary



def collect_data(qc_input):
    ase_atoms = []
    all_data = []
    start, stop = [], []
    energy = []
    for i,  lines in enumerate(qc_input):
        all_data.append(lines.split())
        if "current energy" in lines:
            total = float(lines.split()[2])
            energy.append(total)
        if 'Atoms' in lines:
            start.append(i)
        if 'Lattice vectors (angstrom)' in lines:
            stop.append(i)
    for i, j, energy in zip(start, stop, energy):
        positions = []
        elements = []
        tmp = all_data[i+2:j-1]
        for sys in tmp:
            elements.append(sys[1])
            pos = [float(x) for x in sys[2:]]
            positions.append(pos)
        a = [float(x) for x in all_data[j+1][1:]]
        b = [float(x) for x in all_data[j+2][1:]]
        c = [float(x) for x in all_data[j+3][1:]]
        lattice = [a, b, c]
        ase_atom = Atoms(symbols=elements, positions=positions,
                         cell=lattice, pbc=True)
        ase_atoms.append(ase_atom)
    return ase_atoms, energy


def energy_gradients(data):
    all_gradients = {}
    gradients = data.read_section('History')
    for keys in gradients:
        if 'Gradients' in keys:
            key = keys.replace('Gradients(', '')
            key = int(key.replace(')', ''))
            all_gradients[key] = [gradients[keys][i:i + 3]
                                  for i in range(0, len(gradients[keys]), 3)]

    return all_gradients


def extract(output):
    '''
    '''
    elements, positions, cell = coordlibrary.collect_coords(output)
    New_atoms = Atoms(symbols=elements, positions=positions,
                      cell=cell, pbc=True)
    energy = filetyper.ams_energetic(output)
    return New_atoms, energy


def dict_to_pickle(all_files):
    all_data = {}
    failed = []
    missing = []
    base = os.getcwd()
    for files in all_files:
        tmp = {}
        base_name = files.split('/')[-2]
        output = files+base_name+'.out'
        exp_path = '/'.join(files.split('/')[:-4])
        full_exp_path = exp_path+'/Edited/Valid/'+base_name+'.cif' 
        if os.path.exists(full_exp_path):
            input_structure = read(full_exp_path)
            try:
                new_atoms, energy = extract(output)
                tmp['GFN_atom'] = new_atoms
                tmp['EXP_atom'] = input_structure
                tmp['Energy'] = energy
                all_data[base_name] = tmp
            except:
                failed.append(base_name+'\n')
        else:
            missing.append(base_name+'\n')

    filetyper.put_contents('failed.txt', failed)
    filetyper.put_contents('missing.txt', missing)
    encorder = filetyper.AtomsEncoder
    filetyper.json_to_aseatom(all_data, encorder, 'Optmised_Unoptimised.json')
    return all_data
gfn_files = sorted(glob.glob('/home/mmm0555/Scratch/work/MOF_database/GFN_xTB/MOFData/*/'))