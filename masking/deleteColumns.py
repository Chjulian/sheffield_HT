import sys
import Bio.SeqIO

###############################################################################################
#This script removed columns from an alignment
#It takes 2 arguments 1. the path to the alignment 2. a csv file with the columns to delete in the first column
###############################################################################################


def getErrors(error_file):
    errors = []
    with open(error_file, "r") as infile:
        for line in infile:
            columns = line.split(",")[0]
            if "-" in columns:
                start_col = columns.split("-")[0]
                end_col = columns.split("-")[1]
                errors = errors + [*range(int(start_col), int(end_col) + 1)]
            else:
                errors = errors + [int(columns)]
    return errors

def printNew(aln_file, columns_to_delete):
    in_aln = list(Bio.SeqIO.parse(aln_file, "fasta"))
    for sequence in in_aln:
        as_str = str(sequence.seq)
        print(">" +  sequence.id)
        new_seq = []
        for position in range(0, len(as_str)):
            if position + 1 not in columns_to_delete:
                new_seq.append(as_str[position])
        print("".join(new_seq))

def main():
    start_aln = sys.argv[1]
    aln_errors = sys.argv[2]
    printNew(start_aln, getErrors(aln_errors))

if __name__ == "__main__":
    main()