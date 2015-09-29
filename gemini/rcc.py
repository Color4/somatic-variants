"""
Operations specific to RCC analysis.
"""

import sys
import re
import argparse
import itertools
import pandas as pd

import gemini_operations as gem_ops

def split_id(sample):
    """
    For a full sample name, returns a tuple (sample_name, plate#, tissue type, replicate)
    """
    sample_id, plate_num, tissue, replicate = sample.split('-')
    plate_num = int(plate_num[-1])
    return sample_id, plate_num, tissue, replicate

def get_replicates(gemini_db):
    """
    Returns a dataframe of sample information with the following columns:
    samples, plate_num, tissue, replicate
    """

    samples_df = pd.DataFrame([[sample] + list(split_id(sample)) for sample in gem_ops.get_samples(gemini_db)], \
                              columns=["full_name", "sample", "plate", "tissue", "replicate"])
    return samples_df

def print_variants_sample_pk(dataframe, annotations, samples=None):
    """
    Print variants in table with the sample as primary key, so each row denotes one genotype
    called in one sample.
    """

    if not samples:
        samples = [col[4:] for col in list(dataframe.columns.values) if str(col).startswith("gts.")]

    # Columns to show.
    cols = gem_ops.DEFAULT_VAR_COLS + annotations + ["S_and_N", "S_and_T", "T_and_N"]

    # Print header.
    print "\t".join(["sample", "id", "plate", "tissue", "allele_freq", "alt_depth", "depth"] + cols)

    for sample in samples:
        id, plate, tissue, replicate = split_id(sample)
        variants = dataframe[dataframe["gt_quals." + sample] > 0]
        for i, variant in variants.iterrows():
            print "\t".join([sample, id, str(plate), tissue,
                             "%.3f" % variant["gt_quals." + sample],
                             "%.0f" % variant["gt_alt_depths." + sample],
                             "%.0f" % variant["gt_depths." + sample]
                             ] + list([str(value) for value in variant[cols]])
                            )

def print_variants_variant_pk(dataframe, samples=None):
    """
    Print variants in table with the variant as the primary key, so each row denotes a single
    variant and ALL genotypes called for that variant.
    """

    if samples:
        # Select sample columns.
        # Determine which variants have a genotype in sample_fields
        # Select those variants.
        pass

    # Columns to show.
    cols = gem_ops.DEFAULT_VAR_COLS

    # Get genotype columns and samples list.
    if samples is None:
        gts_col_names = [col for col in list(dataframe.columns.values) if str(col).startswith("gts.")]
        samples = [col[4:] for col in gts_col_names]

    # Print header.
    print "\t".join(cols + ["samples", "allele_freqs", "alt_depths", "depths"])

    # Print variant information.
    for i, variant in dataframe.iterrows():
        # Get genotypes for variant.
        gts = variant[gts_col_names]
        filtered_gts = gts[gts != './.'].index

        # Get allele frequencies, alt depths, and depths.
        gt_quals = filtered_gts.map(lambda n: n.replace("gts", "gt_quals"))
        gt_alt_depths = filtered_gts.map(lambda n: n.replace("gts", "gt_alt_depths"))
        gt_depths = filtered_gts.map(lambda n: n.replace("gts", "gt_depths"))

        # Print variant if there are genotypes called.
        if len(filtered_gts) > 0:
            print '\t'.join( \
                list([str(value) for value in variant[cols]]) + \
                [",".join([gts_col[4:] for gts_col in list(filtered_gts)])] + \
                [",".join(["%.3f" % af for af in variant[gt_quals]])] + \
                [",".join(["%.0f" % ad for ad in variant[gt_alt_depths]])] + \
                [",".join(["%.0f" % d for d in variant[gt_depths]])] )

def add_cols(somatic_vars, samples):
    """
    Add columns to flag variants that are shared between samples.
    """

    # Get column names.
    serum_sample = samples[samples.tissue == 'S'].full_name.item()
    tumor_sample = samples[samples.tissue == 'T'].full_name.item()
    normal_sample = samples[samples.tissue == 'N'].full_name.item()
    ffpe_samples = list(samples[samples.tissue == 'F'].full_name)

    def all_samples_have_gt(row, cols=None):
        """
        Returns true iff all cols denote valid genotypes.
        """
        gt_quals = row[["gt_quals." + col for col in cols]]
        return ( (gt_quals > 0).sum() ) == len(cols)

    # Add columns to flag variants shared between:
    # S and N
    # S and T
    # T and N
    somatic_vars['S_and_N'] = somatic_vars.apply(all_samples_have_gt, axis=1, cols=[serum_sample, normal_sample])
    somatic_vars['S_and_T'] = somatic_vars.apply(all_samples_have_gt, axis=1, cols=[serum_sample, tumor_sample])
    somatic_vars['T_and_N'] = somatic_vars.apply(all_samples_have_gt, axis=1, cols=[tumor_sample, normal_sample])

if __name__ == "__main__":
    # Argument setup and parsing.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("gemini_db", help="Gemini database to use")
    parser.add_argument("--min-allele-freq", type=float, default=0.1, help="Allele frequency to use")
    parser.add_argument("--min-depth", type=int, default=50, help="Minimum depth to use")
    parser.add_argument("--min-alt-depth", type=int, default=5, help="Minimum alternate depth to use")
    args = parser.parse_args()
    gemini_db = args.gemini_db
    min_allele_freq = args.min_allele_freq
    min_depth = args.min_depth
    min_alt_depth = args.min_alt_depth

    annotations = ["TCGA_RCC"]

    # Get genotypes datafame.
    variants = gem_ops.get_genotypes_df(gemini_db, annotations, min_novel_af=min_allele_freq, \
                                        min_anno_af=min_allele_freq, min_alt_depth=min_alt_depth, min_depth=min_depth)

    # Keep only potentially somatic.
    somatic_vars = gem_ops.reduce_to_somatic(variants)

    # Add columns with additional information.
    add_cols(somatic_vars, get_replicates(gemini_db))

    # Print genotypes table.
    print_variants_sample_pk(somatic_vars, annotations)


# This is OLD code not currently used:
"""if __name__ == "__main__":
    gemini_db = sys.argv[1]
    #rvn.print_replicate_tables2(gemini_db, ["TCGA_RCC"], get_replicates)

    # Get somatic variants with replicates in order.
    replicates = sorted( get_replicates(gemini_db, flatten=True) )
    somatic_vars_by_sample = gem_ops.get_somatic_vars_by_sample2(gemini_db, ["TCGA_RCC"], hotspot_af=0.0, samples=replicates)
    somatic_vars_by_sample[["variant_id", "start"]].astype(int)

    for somatic_df in [ somatic_vars_by_sample ]:
        print "********************************* "

        # Print stats.

        print '\t'.join(["Tissue", "Count", "Min", "Max", "Mean", "Median"])
        for tissue in TISSUE_IDS:
            rows_ew = somatic_df["sample"].str.endswith(tissue)
            tissue_rows = somatic_df.loc[rows_ew[rows_ew == True].index]
            stats = ["%.3f" % s for s in [tissue_rows['allele_freq'].min(), tissue_rows['allele_freq'].max(), tissue_rows['allele_freq'].mean(), tissue_rows['allele_freq'].median()] ]
            print '\t'.join([tissue, str(len(tissue_rows))] + stats)

        # Print somatic variants.
        output_file = open("somatic_sample_pk.txt", 'w')
        print somatic_df.to_string(buf=output_file, index=False)
        print somatic_df.to_string(index=False)
        output_file.close()

        # Print gene counts.
        print somatic_df['gene'].value_counts()

        # Print information about each variant: gene, chrom, start, and genotype information.
        output_file = open("somatic_variant_pk.txt", "w")
        temp = sys.stdout
        with output_file as sys.stdout:
            print_variants_table(somatic_df)
        output_file.close()
        sys.stdout = temp
        print_variants_table(somatic_df)

        print ""
        print "Variants shared in replicates"
        replicates = get_replicates(gemini_db)
        for id in sorted( replicates.keys() ):
            reps = replicates[id]
            print ""
            print '**** ', id, ','.join(reps)
            print print_variants_table( somatic_df[somatic_df["sample"].isin(reps)] )
"""
