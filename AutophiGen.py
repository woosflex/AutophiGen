"""
DISCALIMER: This program is part of final assessment for CS50P coursework submitted by:
Name: Adnan Raza
City, Country: New Delhi, India
Program Name: AutophiGen
Description: This program performs basic phylogenetic analysis from searching NCBI databases using BLAST,
aligning BLAST search results using MUSCLE and performing simple phylgenetic analysis. This program is
to be used as a reference generator and is not reliable enough to be included in scientific journals.
(I'm looking forward to working on that aspect of the program in near future)
"""
# Import modules
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ete3 import Tree, TreeStyle
import os
import requests
import sys
from time import sleep
from validator_collection import is_email


# Classes
# Error for format
class FormatError(Exception):
    pass

# Error for BLAST
class BlastError(Exception):
    pass


# Error for MUSCLE
class MuscleError(Exception):
    pass


def main():
    # Create and Check for arguements
    args = apg_args()
    # sequence functionality
    seq_functions(args.input, args.email)


# Main function for arguements
def apg_args():
    args = arg_parser()
    file_check(args.input)
    email_check(args.email)
    return args


# Create arguements
def arg_parser():
    parser = argparse.ArgumentParser(
        description="Automates phylogenetic analysis from Fasta files"
    )
    parser.add_argument(
        "-i", "--input", help="Input file containing sequences", type=str, required=True
    )

    parser.add_argument(
        "-e", "--email", help="Email address for BLAST", type=str, required=True
    )
    return parser.parse_args()


# Check for correct input file
def file_check(file):
    try:
        result = os.path.splitext(file)
        if result[1] not in [".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"]:
            raise FormatError
        return "valid format"
    except FormatError:
        sys.exit("Not a fasta file")


# Check for valid email
def email_check(email):
    try:
        if not is_email(email):
            raise FormatError
        return "valid mail address"
    except FormatError:
        sys.exit("Invalid email")


# main function for BLAST
def seq_functions(file, email):
    # Open file and run through SimpleFastaParser

    with open(file) as input:
        for seq_id, seq in SimpleFastaParser(input):
            try:
                print(f"Sending {seq_id} for BLAST")
                # Make a BLAST query for sequence
                blast(seq, seq_id, email)
                # Align sequences to produce result
                aligned_seq = align(f"result_{seq_id}.fasta", seq_id, email)
                # Performs phylogenetic analysis
                phylo = phylogeny(aligned_seq, seq_id, email)
                tree_viewer(phylo)
            except BlastError:
                print(f"Error in performing BLAST for {seq_id}")
                continue
            except MuscleError:
                print(f"Error in Aligning {seq_id}")


# handles blast functionality
def blast(seq, seq_id, email):
    # Send for BLAST
    NCBIWWW.email = email
    blast_records = NCBIWWW.qblast(
        program="blastn", database="nt", sequence=seq, expect=0.04, hitlist_size=8
    )
    # Parse result and save as Fasta file
    result_parse(seq, blast_records, seq_id)


# parse the result in XML format
def result_parse(seq, result, seq_id):
    # reads the result and parses
    blast_record = NCBIXML.read(result)
    # Create a new file for Blast results
    with open(f"result_{seq_id}.fasta", "w") as output_handle:
        output_handle.write(f">{seq_id}\n{seq}\n")
    # Append the blast result in newly created file
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            with open(f"result_{seq_id}.fasta", "a+") as output:
                if "complete" in alignment.title:
                    continue
                output.write(f">{alignment.title}\n{hsp.sbjct}\n")


# align the result file of BLAST
def align(seqs, seq_id, email):
    print(f"Aligning BLAST results of {seq_id}")
    # Open BLAST result and store in string format
    sequence_str = ""
    with open(seqs) as sequences:
        sequences = SeqIO.parse(sequences, "fasta")
        seq_count = 0
        for sequence in sequences:
            sequence_str += f">{sequence.description}\n{sequence.seq}\n"
            seq_count += 1
        if seq_count == 1:
            raise BlastError
    # Open result file and send sequences to Muscle API
    send = send_sequences(sequence_str, email, seq_id)
    print(f"Job ID: {send.content.decode()}")
    # Check if process is finished
    check_status(send.content)
    # Recieve results
    result = get_result(send.content)
    return result.content.decode()


# Send sequences for alignment to MUSCLE API
def send_sequences(seqs, email, seq_id):
        headers = {
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "text/plain",
        }

        data = {
            "email": email,
            "sequence": seqs,
            "title": seq_id,
        }

        send = requests.post(
            "https://www.ebi.ac.uk/Tools/services/rest/muscle/run",
            headers=headers,
            data=data,
        )
        if str(send) != "<Response [200]>":
            print("Response not 200")
            raise MuscleError
        
        print(f"Aligning {seq_id}")

        return send


# Check status of MUSCLE query with interval of 5 sec
def check_status(job_id):
    while True:
        s_header = {"Accept": "text/plain"}
        status = requests.get(
            f"https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{job_id.decode()}",
            headers=s_header,
        )
        print(status.content.decode())
        if "FINISHED" in status.content.decode():
            break
        sleep(5)


# Get result of query from SimplePhylo API
def get_result(job_id):
    r_headers = {"Accept": "text/plain"}
    return requests.get(
        f"https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{job_id.decode()}/fa",
        headers=r_headers,
    )


# phylogenetic analysis
def phylogeny(aligned_seqs, seq_id, email):
    print(f"Performing phylogenetic analysis on {seq_id}")
    send = send_sequences_phylo(aligned_seqs, email, seq_id)
    print(f"Job ID: {send.content.decode()}")
    check_status_phylo(send.content.decode())
    result = get_result_phylo(send.content.decode())
    return result.content.decode()


# send sequences for phylogenetic analysis on EMBL-EBI API
def send_sequences_phylo(aligned_seqs, email, seq_id):
    headers = {
        "Content-Type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
    }

    data = {
        "email": email,
        "sequence": aligned_seqs,
        "title": seq_id,
    }

    send = requests.post(
        "https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/run",
        headers=headers,
        data=data,
    )

    return send


# Check status of Phylo query with interval of 5 sec
def check_status_phylo(job_id):
    while True:
        sr_header = {"Accept": "text/plain"}
        status = requests.get(
            f"https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/status/{job_id}",
            headers=sr_header,
        )
        print(status.content.decode())
        if "FINISHED" in status.content.decode():
            break
        sleep(5)


# Get result of query from MUSCLE API
def get_result_phylo(job_id):
    pr_headers = {"Accept": "application/json"}
    return requests.get(
        f"https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/result/{job_id}/tree",
        headers=pr_headers,
    )


# View phylogenetic tree
def tree_viewer(tree_data):
    gene_tree_nw = Tree(tree_data)
    ts = TreeStyle()
    ts.branch_vertical_margin = 10
    gene_tree_nw.show(tree_style=ts)


if __name__ == "__main__":
    main()
