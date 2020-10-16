with open('GCF_000859885.1_ViralProj15197_protein.faa') as f:
   text = f.read()
   for protein in text.split('>'):
     if not protein:
       continue
     with open('smallpox/' + protein[:11], 'w+') as f2:
       f2.write('>' + protein)