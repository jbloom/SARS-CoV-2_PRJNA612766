"""Generate pileup and nucleotide identity calls from BAM file.

For command line usage, see::

    python bam_pileup.py --help

"""

__author__ = 'Jesse Bloom'
__email__ = 'jbloom@fredhutch.org'


import argparse
import os

import Bio.SeqIO

import pandas as pd

import pysam


def bam_pileup(bam,
               ref,
               ref_fasta,
               pileup_csv,
               minq=20,
               bai=None,
               add_cols=None,
               ):
    """Pileup and nucleotide identity calls from BAM."""
    if not os.path.isfile(bam):
        raise IOError(f"cannot find `bam` {bam}")
    if bai is None:
        bai = bam + '.bai'
    if not os.path.isfile(bai):
        raise IOError(f"cannot find `bai` {bai}")

    with pysam.AlignmentFile(bam, 'rb', index_filename=bai) as bamfile:
        if ref not in bamfile.references:
            raise ValueError(f"`ref` {ref} not among valid references:\n\t" +
                             '\n\t'.join(bamfile.references))
        count_df = pd.DataFrame({nt: counts for nt, counts in
                                 zip('ACGT',
                                     bamfile.count_coverage(
                                            contig=ref,
                                            quality_threshold=minq)
                                     )
                                 })
        if len(count_df) != bamfile.get_reference_length(ref):
            raise ValueError('not expected number of sites')

    ref_seq = str(Bio.SeqIO.read(ref_fasta, 'fasta').seq)
    if len(ref_seq) != len(count_df):
        raise ValueError('sequence in `ref_fasta` not correct length')

    nt_cols = ['A', 'C', 'G', 'T']
    assert count_df.columns.tolist() == nt_cols
    count_df = (
        count_df
        .assign(depth=lambda x: x[nt_cols].sum(axis=1),
                site=lambda x: x.index + 1,
                ref_nt=lambda x: list(ref_seq.upper()),
                mut_depth=lambda x: x.apply(
                                lambda r: sum(r[nt] for nt in nt_cols
                                              if nt != r['ref_nt']),
                                axis=1)
                )
        [['site', 'ref_nt', 'depth', 'mut_depth', *nt_cols]]
        )
    for col_name, col_val in add_cols:
        if col_name in set(count_df.columns):
            raise ValueError(f"column {col_name} already in output CSV")
        count_df[col_name] = col_val

    count_df.to_csv(pileup_csv, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                description='Pileup and nucleotide identity calls from BAM.',
                )
    parser.add_argument('--bam',
                        required=True,
                        help='BAM file',
                        )
    parser.add_argument('--ref',
                        required=True,
                        help='name of contig over which we get coverage',
                        )
    parser.add_argument('--ref_fasta',
                        required=True,
                        help='FASTA file with reference',
                        )
    parser.add_argument('--pileup_csv',
                        required=True,
                        help='name of created pileup CSV',
                        )
    parser.add_argument('--minq',
                        help='only count bases with Q >= this',
                        )
    parser.add_argument('--bai',
                        help='BAI file, if none then BAM file suffixed by .bai'
                        )
    parser.add_argument('--add_cols',
                        nargs=2,
                        metavar=('col_name', 'col_value'),
                        action='append',
                        help='columns to add to pileup CSV, can use >1 times',
                        )
    args = vars(parser.parse_args())
    bam_pileup(**args)

