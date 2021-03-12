/usr/bin/time -a -o howlongittook.txt \
	python build_db_from_frags.py --cloud \
	--known-seq-path ~/blastdb/known-seqs-db/ \
	resources/hazard_db/hazards.fragset.json \
	> trace_build_db.txt
