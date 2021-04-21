import argparse

import Bio.SeqIO


def strip_fasta_to_head_id(fasta):
    """Strip FASTA header to just retain first space-delimted word."""
    seqs = list(Bio.SeqIO.parse(fasta, 'fasta'))
    with open(fasta, 'w') as f:
        for s in seqs:
            f.write(f">{s.id.split()[0]}\n{str(s.seq)}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='strip FASTA header to id (first space-delimited word)')
    parser.add_argument('--fasta',
                        required=True,
                        help='name of FASTA file')
    args = vars(parser.parse_args())
    strip_fasta_to_head_id(**args)
