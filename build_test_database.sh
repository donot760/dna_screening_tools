python funtrp_sequence_struct_pick.py &&
echo done picking fragments &&
python seq_file_to_frag_set.py resources/known_seqs/e_coli/e_coli_genbank.faa resources/known_seqs/e_coli_aa.fragset.json &&
echo done splitting e. coli into fragments &&
python seq_file_to_frag_set.py resources/known_seqs/taterapox/taterapox_aas_and_dna_genbank.gb resources/known_seqs/taterapox_aa.fragset.json &&
echo done splitting taterapox into fragments &&
python seq_file_to_frag_set.py resources/known_seqs/rinderpest/rinderpest_aas_and_dna.gb resources/known_seqs/rinderpest_aa.fragset.json &&
echo done splitting rinderpest into fragments &&
python build_db_from_frags.py resources/hazard_db/hazards.fragset.json &&
echo done building database &&
echo finished all steps successfully.
