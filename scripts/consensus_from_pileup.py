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
    cols = ['site', 'depth'] + nts
    pileup = pd.read_csv(pileup)
    if not set(cols).issubset(set(pileup.columns)):
        raise ValueError(f"{pileup} lacks columns {nts}")
    if len(pileup) != len(pileup[cols].drop_duplicates()):
        raise ValueError(f"duplicated columns in {pileup}")

    first_site = pileup['site'].min()
    last_site = pileup['site'].max()
    if set(pileup['site']) != set(range(first_site, last_site + 1)):
        raise ValueError(f"{pileup} does not contain sequential sites")

    site_to_nt = (
        pileup
        .melt(id_vars=['site', 'depth'],
              value_vars=nts,
              var_name='nt',
              value_name='count',
              )
        .assign(frac=lambda x: x['count'] / x['depth'])
        .query('count > @min_coverage')
        .query('frac > @min_frac')
        .sort_values(['site', 'count'], ascending=[True, False])
        .groupby('site')
        .aggregate({'nt': 'first'})
        ['nt']
        .to_dict()
        )

    seq = ''.join([site_to_nt[r] if r in site_to_nt else 'N'
                   for r in range(first_site, last_site + 1)])
    with open(consensus, 'w') as f:
        f.write(f">{fasta_header}\n{seq}")


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
