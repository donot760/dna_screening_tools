import sys, os, re

# from django -> utils -> text.py
def valid_filename(s):
    """
    Return the given string converted to a string that can be used for a clean
    filename. Remove leading and trailing spaces; convert other spaces to
    underscores; and remove anything that is not an alphanumeric, dash,
    underscore, or dot.
    >>> get_valid_filename("john's portrait in 2004.jpg")
    'johns_portrait_in_2004.jpg'
    """
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'(?u)[^-\w.]', '', s)

_, in_file, out_path = sys.argv
with open(in_file) as f:
   text = f.read()
   for protein in text.split('>'):
     if not protein:
       continue
     with open(os.path.join(out_path, valid_filename(protein.split()[0] +
             '_' + os.path.basename(in_file).split('.')[0]).replace('.', '_') +
             '_split.faa'), 'w+') as f2:
       f2.write('>' + protein)
