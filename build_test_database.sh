python funtrp_sequence_pick.py --noplot &&
python seq_file_to_frag_set.py resources/genbank/e_coli/e_coli_genbank.faa resources/genbank/e_coli_aa.fragset.json &&
python build_db_from_frags.py resources/hazard_db/hazards.fragset.json
