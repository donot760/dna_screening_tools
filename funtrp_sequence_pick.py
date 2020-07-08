# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:44:50 2020

@author: jave3
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def list_files(directory):
#    os.chdir(directory)
    current_dir = os.getcwd()

    new_dir = "resources/select_agent_hazard/funtrp_results/"
    os.chdir(new_dir)
    csv_files = [fname for fname in os.listdir() if fname.split('.')[-1] == 'csv'] # gets all the files in new_dir with a .csv extension
    os.chdir(current_dir) # back to original directory to access output text file without convoluted relative paths
    return new_dir, csv_files

def residue_and_toggles(filename):
    untouched_toggles = []
    residues = []
    rheos = []
    neutrals = []
    
    with open(filename) as csv_file:
        csv_reader  = csv.reader(csv_file, delimiter = ',')
        next(csv_reader)#skip the header
        for line in csv_reader:#this goes through every line
            prediction_t = float(line[5])#prediction toggle
            prediction_r = float(line[9])
            prediction_n = float(line[6])
            residue = line[11]
            untouched_toggles.append(prediction_t)#these prediction toggles are straight from the data
            residues.append(residue)
            rheos.append(prediction_r)
            neutrals.append(prediction_n)

    return residues, untouched_toggles, rheos, neutrals

def shift_toggles(toggle_values, shift):
    
    shifted_toggles = list(map(lambda toggle: toggle + shift, toggle_values))
#    print(shifted_toggles)
    return shifted_toggles

def create_possible_sequences(residues, toggles, rheos, neutrals, window = 19):
    count = 0
    shifted_toggles = shift_toggles(toggles, .3)#changed to have the shift here, enter origianl toggle values instead
    possible_tval_sequences = []
    possible_residue_sequences = []
    t_r_n = []#toggles, neutrals, and rheos
    for toggle in shifted_toggles:#error: going to have the same number in the data, will think that something with an index of 130 is actually a value with index of 2
        if (count + window) <= len(shifted_toggles):
            possible_tval_sequences.append(shifted_toggles[count: count + window])
            possible_residue_sequences.append(residues[count: count + window])
            t_r_n.append(toggles[count: count + window] + rheos[count: count + window] + neutrals[count: count + window])
            
            
        count += 1
    
    for index in range(len(possible_tval_sequences)):
        possible_tval_sequences[index] = np.prod(possible_tval_sequences[index])#get the product for each sequence
        
    max_vals = []
    
    copied_sequences = possible_tval_sequences.copy()
    for index in range(10):
        max_ = 0
        for i in range(len (copied_sequences)):
            if copied_sequences[i] > max_:
                max_ = copied_sequences[i]
                
        copied_sequences.remove(max_)
        max_vals.append(max_)
            
    print(max_vals)
    index_values = []
    for value in max_vals:
        index_values.append(possible_tval_sequences.index(value))
            
#    print(index_values)
#    
#    for i in index_values:
#        if possible_residue_sequences
    filtered_sequence = []
    for residue_sequence in possible_residue_sequences:
        if possible_residue_sequences.index(residue_sequence) in index_values:
            filtered_sequence.append(residue_sequence)
#            possible_residue_sequences.remove(residue_sequence)
#    print(filtered_sequence)
    
    string_list = []
    with open("resources/aa_fragment_picks.txt","w+") as f:
        for sequence in filtered_sequence:
            f.write("{}\r\n".format(sequence))
    
    
    with open("resources/aa_fragment_picks.txt","r") as f:
        lines = f.readlines()
#        print(lines)
    
    
    for line in lines:
        new_string = ""
        for aa in line:
            if aa.isalpha():
                new_string += aa
        string_list.append(new_string)
#    print(string_list)

    with open("resources/aa_fragment_picks.txt","w+") as f:
        for string in string_list:
            f.write("{}\r\n".format(string))
   
        
    return possible_tval_sequences, possible_residue_sequences

#def products_for_sequences(list_of_sequences, residues):
#    
#    for index in range(len(list_of_sequences)):
#        list_of_sequences[index] = np.prod(list_of_sequences[index])#get the product for each sequence
#        
#    return list_of_sequences

#def graph(sequence_products):
#    
#    plt.plot(range(len(sequence_products)), sequence_products)
#    plt.show()
    

    
if __name__ == '__main__':
    directory = list_files("UROP 2020")
#    files = os.listdir(directory)[1:]
#    new_dir = os.chdir(os.getcwd() + "/../funtrp_results/")
#    os.getcwd()
#    list_files(new_dir)
    files = directory[1]
    for file in files:

#        print("current csv: ", end = '')
#        print(file)
        r_t = residue_and_toggles(directory[0] +"/"+ str(file))
#    r_t = residue_and_toggles(os.getcwd() + "/../funtrp_results/NP_042046_funtrp.csv")
#    shifted_values = (shift_toggles(r_t[1], .3))
        tval_possibilities = create_possible_sequences(r_t[0], r_t[1] ,r_t[2], r_t[3])
#    products = products_for_sequences(tval_possibilities[0], tval_possibilities[1])
    
#    graph(products)
