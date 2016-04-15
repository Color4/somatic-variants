#
# Simple script to count shared and unique variants between different tissues.
# Input file was generated using:
#   python rcc.py augment_vars --results-file find_vars_results_all_minaf0.02_ad5_d20.txt --tissue N,S
# Usage of this script:
#   python normal_vs_serum.py <input_tsv>
#

import sys
import itertools
import pandas as pd


def samples_shared(variants_df):
    """
    Compare samples from same individual pairwise and count shared variants.
    """

    # Iterate through each individual.
    for name, id_group in variants_df.groupby(["id"]):
        print name, "*********"

        # Compare each combination of samples.
        for collection1, collection2 in itertools.combinations(id_group.groupby(["tissue", "replicate"]), 2):
            (tissue1, rep1), sample_variants1 = collection1
            (tissue2, rep2), sample_variants2 = collection2
            print "\t%s%s: %i variants, %s%s: %i variants" % (tissue1, rep1, sample_variants1["id"].count(), tissue2, rep2, sample_variants2["id"].count())
            shared = sample_variants1.merge(sample_variants2, on=["variant_id"])
            print "\t\tShared: ", len(shared)

            # If comparing serum and normal, output unique variants.
            # if (tissue1 == "S" and tissue2 == "N"):
            #     print "@@@@@@", len(sample_variants1[~sample_variants1["variant_id"].isin(shared["variant_id"])])

            if (tissue1 == "N" and tissue2 == "S"):
                output = open("%s_%s%s_unique_vs_%s%s.txt" % (name, tissue2, rep2, tissue1, rep1), "w")
                # print "$$$$$$", len(sample_variants2[~sample_variants2["variant_id"].isin(shared["variant_id"])])
                serum_unique = sample_variants2[~sample_variants2["variant_id"].isin(shared["variant_id"])]
                output.write( serum_unique.to_csv(sep="\t") )
                output.close()


if __name__ == "__main__":
    infile = sys.argv[1]

    variants_df = pd.read_csv(infile, sep="\t")
    variants_df = variants_df[ (variants_df["allele_freq"] >= 0.1) & (variants_df["alt_depth"] >= 10)]
    samples_shared(variants_df)
