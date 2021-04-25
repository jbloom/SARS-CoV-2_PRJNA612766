"""Call consensus sequence from pileup CSV."""

import argparse

import pandas as pd


def consensus_from_pileup(pileup,
                          consensus,
                          fasta_header,
                          min_coverage,
                          min_frac,
                          ):
    """Call consensus from pileup CSV."""
    nts = ['A', 'C', 'G', 'T']
    cols = ['site', *nts]
    pileup = pd.read_csv(pileup)
    if not set(cols).issubset(set(pileup.columns)):
        raise ValueError(f"{pileup} lacks columns {nts}")
    if len(pileup) != len(pileup[cols].drop_duplicates()):
        raise ValueError(f"duplicated rows in {pileup}")

    first_site = pileup['site'].min()
    last_site = pileup['site'].max()
    if set(pileup['site']) != set(range(first_site, last_site + 1)):
        raise ValueError(f"{pileup} does not contain sequential sites")

    seq = ''.join(
        pileup
        .sort_values('site')
        .assign(depth=lambda x: x[nts].sum(axis=1),
                most_common=lambda x: x[nts].idxmax(axis=1),
                most_common_depth=lambda x: x[nts].max(axis=1),
                most_common_frac=lambda x: x['most_common_depth'] / x['depth'],
                meets_criteria=lambda x: (
                            (x['most_common_depth'] >= min_coverage) &
                            (x['most_common_frac'] > min_frac)),
                consensus=lambda x: x['most_common'].where(x['meets_criteria'],
                                                           'N')
                )
        ['consensus']
        .tolist()
        )

    with open(consensus, 'w') as f:
        f.write(f">{fasta_header}\n{seq}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Call consensus from pileup CSV.')
    parser.add_argument('--pileup',
                        required=True,
                        help='input pileup CSV file')
    parser.add_argument('--consensus',
                        required=True,
                        help='created FASTA with consensus')
    parser.add_argument('--fasta_header',
                        required=True,
                        help='header for FASTA')
    parser.add_argument('--min_coverage',
                        type=int,
                        required=True,
                        help='require >= this coverage to call a site')
    parser.add_argument('--min_frac',
                        type=float,
                        required=True,
                        help='require most common base > this fraction')
    args = vars(parser.parse_args())
    consensus_from_pileup(**args)
