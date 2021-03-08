# Command to download the whole fasta for the db with name given by the first argument:
echo "Downloading all blast database files for database named $1. Completion time:"
time docker run --rm \
     -v $HOME/blastdb:/blast/blastdb:rw \
     -w /blast/blastdb \
     ncbi/blast \
     update_blastdb.pl "$1"
echo "Done."
echo "Converting database named $1 to a raw fasta file $1.fsa. Completion time:"
time docker run --rm \
     -v $HOME/blastdb:/blast/blastdb:rw \
     -w /blast/blastdb \
     ncbi/blast \
     blastdbcmd -entry all -db "$1" \
     -out "$1.fsa"
echo "Done."

