# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:44:50 2020

@author: jave3
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

def list_files(directory):
    pass

def residue_and_toggles(filename):
    untouched_toggles = []
    residues = []
    
    with open(filename) as csv_file:
        csv_reader  = csv.reader(csv_file, delimiter = ',')
        next(csv_reader)#skip the header
        for line in csv_reader:#this goes through every line
            prediction_t = float(line[5])#prediction toggle
            residue = line[11]
            untouched_toggles.append(prediction_t)#these prediction toggles are straight from the data
            residues.append(residue)
#    print(residues, untouched_toggles)
    return residues, untouched_toggles

def shift_toggles(toggle_values, shift):
    
    shifted_toggles = list(map(lambda toggle: toggle + 0.3, toggle_values))
#    print(shifted_toggles)
    return shifted_toggles

def create_possible_sequences(shifted_toggles, residues, window = 19):
    count = 0
    possible_tval_sequences = []
    possible_residue_sequences = []
    for toggle in shifted_toggles:#error: going to have the same number in the data, will think that something with an index of 130 is actually a value with index of 2
        if (count + window) <= len(shifted_toggles):
            possible_tval_sequences.append(shifted_toggles[count: count + window])
            possible_residue_sequences.append(residues[count: count + window])
        count += 1
            
#    string_list = []
#    with open("resources/aa_fragment_picks.txt","w+") as f:
#        for sequence in possible_residue_sequences:
#            f.write("{}\r\n".format(sequence))
#    
#    
#    with open("resources/aa_fragment_picks.txt","r") as f:
#        lines = f.readlines()
#        print(lines)
#    
#    
#    for line in lines:
#        new_string = ""
#        for aa in line:
#            if aa.isalpha():
#                new_string += aa
#        string_list.append(new_string)
##    print(string_list)
#
#    with open("resources/aa_fragment_picks.txt","w+") as f:
#        for string in string_list:
#            f.write("{}\r\n".format(string))
   
        
    return possible_tval_sequences, possible_residue_sequences

def products_for_sequences(list_of_sequences):
    
    for index in range(len(list_of_sequences)):
        list_of_sequences[index] = np.prod(list_of_sequences[index])#get the product for each sequence
        
    string_list = []
    with open("resources/aa_fragment_picks.txt","w+") as f:
        for sequence in list_of_sequences:
            f.write("{}\r\n".format(sequence))
    
    
    with open("resources/aa_fragment_picks.txt","r") as f:
        lines = f.readlines()
        print(lines)
    
    
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
    
    return list_of_sequences

def graph(sequence_products):
    
    plt.plot(range(len(sequence_products)), sequence_products)
    plt.show()
    
#def filter_(threshold, final_products):
#    copied = final_products.copy()
#    count = 0
#    while count < len(copied):
#        if final_products[count] < threshold:
#            final_products.remove(final_products[count])
#            count += 1
#    
#    print(final_products)
            
    
if __name__ == '__main__':
    r_t = residue_and_toggles("resources/select_agent_hazard/funtrp_results/NP_042046_funtrp.csv")
    shifted_values = (shift_toggles(r_t[1], .3))
    tval_possibilities = create_possible_sequences(shifted_values, r_t[0])
    products = products_for_sequences(tval_possibilities[0])

    graph(products)
