from Bio.Seq import Seq
import os
import time

#from Bio.Alphabet import ExtendedIUPACProtein

genbank_path = os.path.join('resources', 'genbank')


def get_all_files():
    #current_directory = os.getcwd()
    #inside_directory = os.listdir()
    files = []
    for file_ in os.listdir(genbank_path):
        if file_.split(".")[-1] == "fna":
            files.append(os.path.join(genbank_path, file_))
    return files[:3]

def get_all_files_gen():
    files = []
    for file_ in os.listdir(genbank_path):
        if file_.split(".")[-1] == "fna":
            yield os.path.join(genbank_path, file_)
    
def multiple_of_3(reading_frame):
    remainder = len(reading_frame)%3
    
    return reading_frame[:len(reading_frame)-remainder]
    


def get_dna(fname):     
    with open(fname, "r") as f:
        line_list = f.readlines()[1:]
        for line in line_list:
            if line.startswith(">"):#this does not work if you have consecutive lines starting with ">"
                line_list.remove(line)
        whole_sequence = "".join(line_list)
        whole_sequence = whole_sequence.replace('\n', '')
    
    #reading_frame_1 = multiple_of_3(Seq(whole_sequence))
    ##reading_frame_2 = multiple_of_3(Seq(whole_sequence[1:]))
    #reading_frame_3 = multiple_of_3(Seq(whole_sequence[2:]))

    
    return Seq(whole_sequence)
    
def reading_frame(whole_seq):
    reading_frame_1 = multiple_of_3(whole_seq)
    reading_frame_2 = multiple_of_3(whole_seq[1:])
    reading_frame_3 = multiple_of_3(whole_seq[2:])
    
    return reading_frame_1, reading_frame_2, reading_frame_3 



def make_codons(sequence):
    nucleotide_list = list(sequence)
    count = 0
    codons = []
    remainder = len(nucleotide_list) % 3
    
    while True: # for nucleotide in nucleotide_list[:len(nucleotide_list) - remainder]:#this was used in case % 3 != 0
    
        try: 
            codon = nucleotide_list[count] + nucleotide_list[count+1] + nucleotide_list[count+2]
            count += 3
            codons.append(codon)
            
        except IndexError:
            break
            
    return codons ##returns list of dna sequence separated into codons
    
def make_amino_acids(codons):
    amino_acids = []
    for codon in codons:
        amino_acid = Seq(codon).translate()
        amino_acids.append(amino_acid)
    
    return amino_acids
    
def aa_windows_generator(amino_acids, window = 19):
    windows = []
    for index in range(len(amino_acids)):
        if index + window <= len(amino_acids):
            #element = windows.append(amino_acids[index:index+window])
            yield amino_acids[index:index+window]
            
    
def aa_windows(amino_acids, window = 19):
    windows = []
    for index in range(len(amino_acids)):
        if index + window <= len(amino_acids):
            windows.append(amino_acids[index:index+window])
            
    return windows

def make_protein_seq(amino_acids): #aa is for now, should be codons
  
    methionine_indexes = []#list that holds the indexes for every methionine
    stop_indexes = []#holds indexes for every STOP
    for index in range(len(amino_acisds)):
        if amino_acids[index] == "M":
            methionine_indexes.append(index)
        if amino_acids[index] == "*":
            stop_indexes.append(index)
    start_and_stop_indexes = []#list of tuples of the start and stop indexes
    for index in range(len(methionine_indexes)):#my reasoning is that the num of start codons = number of stops
        start_and_stop_indexes.append((methionine_indexes[index], stop_indexes[index]))
        #my reasoning is that the first start goes with the first stop, second start goes with the second stop, might not be the case
    proteins = []
    for index_tuple in start_and_stop_indexes:
        start = index_tuple[0]
        stop = index_tuple[1]
        proteins.append(amino_acids[start:stop])

    return proteins

if __name__ == '__main__':
    
    #################now does it all work with the file??
    #print("hi")
    #sequence = get_dna("resources/GCF_000008865.2_ASM886v2_genomic.fna")
    #trnsltd_seq_1 = sequence[0].translate()
    
    #seq_1 = "FVSSAGYSSTVFYGDRKVT"   
    #seq_2 = "DPCLSPCTKLKSKWIKDLH"
    
    #first_frame = aa_windows_generator(trnsltd_seq_1) 
    #for window in first_frame:
    #    print(window)
    
    
    #whole_thing = ''
    #for reading_frame in sequence:
    #    translated_frame = reading_frame.translate()
    #    whole_thing += translated_frame
    #windows = aa_windows(whole_thing)
    #window_set = set(windows)
    ###############################################doing it with the set
    
    for fna in get_all_files_gen():
        
     
        print(fna)
        combined_frame = ''
        sequence = get_dna(fna)
        frames = reading_frame(sequence)
        print("in the file")
        translated_1 = frames[0].translate()
        print('frame 1: ' + translated_1[:10])
        translated_2 = frames[1].translate()
        print('frame_2: ' + translated_2[:10])
        translated_3 = frames[2].translate()
        print('frame_3: ' + translated_3[:10])
        combined_frame += translated_1 + translated_2 + translated_3
         
        
        windows = aa_windows_generator(combined_frame)
        #print(next(windows))
        #print(next(windows))
        #print(next(windows))
        for ele in windows:
            print(ele)
       
        #print(seq_1 in window_set)
        #print(seq_2 in window_set)
        
    #####################################################
    

    
            
   
  
        