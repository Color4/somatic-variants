"""
Operations specific to RCC analysis.
"""

import sys
import re
import time
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

def print_variants_sample_pk(dataframe, annotations, joint_cols=None, samples=None):
    """
    Print variants in table with the sample as primary key, so each row denotes one genotype
    called in one sample.
    """

    if not samples:
        samples = [col[4:] for col in list(dataframe.columns.values) if str(col).startswith("gts.")]

    # Columns to show.
    cols = gem_ops.DEFAULT_VAR_COLS + annotations + joint_cols

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

def add_joint_cols(somatic_vars):
    """
    Add columns to dataframe flag variants that are shared/joint between samples. Returns
    names of columns added to dataframe.
    """

    def joint_genotypes(variant_df, all_have_gt_samples=None, any_has_gt_samples=None, min_count=1):
        """
        Returns true for variant iff (a) all samples in all_have_gt_cols have a genotype and (b) at least min_count
        samples in any_has_gt has genotype.
        """

        # Check (a) criterion.
        all_have_gt = True
        if all_have_gt_samples != None:
            all_have_gt = len(variant_df[variant_df["sample"].isin(all_have_gt_samples)]) == len(all_have_gt_samples)

        # Check (b) criterion.
        any_has_gt = True
        if any_has_gt_samples != None:
            any_has_gt = len(variant_df[variant_df["sample"].isin(any_has_gt_samples)]) >= min_count

        return all_have_gt and any_has_gt

    # Add columns for joint variants.
    # NOTE: should set max serum/tumor/normal by looking at data frame and grouping
    # by ID + tissue.
    num_serum_samples = 2
    num_tumor_samples = 2
    num_normal_samples = 2
    for s in range(num_serum_samples):
        for t in range(num_tumor_samples):
            for n in range(num_normal_samples):
                sn_col = 'j_S%i_N%i' % (s, n)
                st_col = 'j_S%i_T%i' % (s, t)
                tn_col = 'j_T%i_N%i' % (t, n)
                stn_col = 'j_S%iT%iN%i' % (s, t, n)
                sf_col = 'j_S%iF' % (s)
                stnf_col = 'j_S%iT%iN%iF' % (s, t, n)
                somatic_vars[sn_col] = False
                somatic_vars[st_col] = False
                somatic_vars[tn_col] = False
                somatic_vars[stn_col] = False
                somatic_vars[sf_col] = False
                somatic_vars[stnf_col] = False

    for an_id, id_vars in somatic_vars.groupby("id"):
        # Get samples.
        serum_samples = list(id_vars[id_vars.tissue == 'S']['sample'].unique())
        tumor_samples = list(id_vars[id_vars.tissue == 'T']['sample'].unique())
        normal_samples = list(id_vars[id_vars.tissue == 'N']['sample'].unique())
        ffpe_samples = list(id_vars[id_vars.tissue == 'F']['sample'].unique())

        # Add columns to flag variants shared between samples.
        # FIXME: there must be a function to replace these loops.
        columns_added = []
        for variant_id, var_group in id_vars.groupby("variant_id"):
            for n, normal_sample in enumerate(normal_samples):
                for s, serum_sample in enumerate(serum_samples):
                    # S and N
                    sn_col = 'j_S%i_N%i' % (s, n)
                    somatic_vars.loc[var_group.index, sn_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, normal_sample])

                    # S and F
                    sf_col = 'j_S%iF' % (s)
                    somatic_vars.loc[var_group.index, sf_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample], any_has_gt_samples=ffpe_samples)

                    for t, tumor_sample in enumerate(tumor_samples):
                        # T and N
                        # S and T
                        # S, T, and N
                        # S, T, N, and F
                        tn_col = 'j_T%i_N%i' % (t, n)
                        st_col = 'j_S%i_T%i' % (s, t)
                        stn_col = 'j_S%iT%iN%i' % (s, t, n)
                        stnf_col = 'j_S%iT%iN%iF' % (s, t, n)
                        somatic_vars.loc[var_group.index, st_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, tumor_sample])
                        somatic_vars.loc[var_group.index, tn_col] = joint_genotypes(var_group, all_have_gt_samples=[tumor_sample, normal_sample])
                        somatic_vars.loc[var_group.index, stn_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, normal_sample, tumor_sample])
                        somatic_vars.loc[var_group.index, stnf_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, tumor_sample, normal_sample], any_has_gt_samples=ffpe_samples)

                if len(serum_samples) == 0:
                    for t, tumor_sample in enumerate(tumor_samples):
                        tn_col = 'j_T%i_N%i' % (t, n)
                        somatic_vars.loc[var_group.index, tn_col] = joint_genotypes(var_group, all_have_gt_samples=[tumor_sample, normal_sample])

    return somatic_vars

if __name__ == "__main__":
    # Argument setup and parsing.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("operation", help="Operation to perform")
    parser.add_argument("--gemini-db", help="Gemini database to use")
    parser.add_argument("--min-allele-freq", type=float, default=0.05, help="Allele frequency to use")
    parser.add_argument("--min-depth", type=int, default=20, help="Minimum depth to use")
    parser.add_argument("--min-alt-depth", type=int, default=5, help="Minimum alternate depth to use")
    parser.add_argument("--sample-pattern", default=".*", help="Regular expression to use to select samples")
    parser.add_argument("--annotations", default="TCGA_RCC", help="Annotations to use as hotspots")
    parser.add_argument("--results-file", default="", help="File of raw somatic results")
    args = parser.parse_args()

    gemini_db = args.gemini_db
    operation = args.operation
    min_allele_freq = args.min_allele_freq
    min_depth = args.min_depth
    min_alt_depth = args.min_alt_depth
    sample_pattern = args.sample_pattern
    annotations = args.annotations.split(",")
    results_file = args.results_file

    if operation == "find_somatic":
        # Output somatic variants.
        samples = [sample for sample in gem_ops.get_samples(gemini_db) if re.search(sample_pattern, sample) > 0]
        print >> sys.stderr, "Samples to process: ", ", ".join(samples)
        all_somatic = []
        for sample in samples:
            start = time.time()
            somatic_vars = gem_ops.get_somatic_vars_in_sample(gemini_db, annotations, sample, \
                                                              min_anno_af=min_allele_freq, min_novel_af=min_allele_freq, \
                                                              min_alt_depth=min_alt_depth, min_depth=min_depth)
            end = time.time()
            all_somatic.append(somatic_vars)
            print >> sys.stderr, sample, len(somatic_vars), end-start

        # Combine all somatic variants together and add id, plate, tissue.
        all_somatic_df = pd.concat(all_somatic)
        all_somatic_df.reset_index(inplace=True, drop=True)

        # Using sample column, create and populate columns for id, plate, tissue, and replicate.
        sample_attrs = pd.DataFrame(list(all_somatic_df["sample"].apply(lambda s: split_id(s))), columns=["id", "plate", "tissue", "replicate"])
        for i, col in enumerate(["id", "plate", "tissue", "replicate"]):
            all_somatic_df.insert(i+1, col, pd.Series())
            all_somatic_df[col] = sample_attrs[col]

        print all_somatic_df.to_csv(sep="\t", index=False, float_format='%.3f'),

    elif operation == "augment_somatic":
        # Read results into dataframe.
        results_df = gem_ops.convert_cols( pd.read_csv(results_file, sep="\t") )

        # Update num_het, add num_het_by_id
        results_df = gem_ops.update_num_het(results_df)
        results_df = gem_ops.update_num_het_by_id(results_df)

        # Add joint columns.
        results_df = add_joint_cols(results_df)

        # Output sample statistics.
        #print results_df.groupby(["sample"]).size()

        #vhl_muts = results_df[(results_df["gene"] == "VHL") & (results_df["alt_depth"] > 10) & (results_df["depth"] > 40)]
        # for tissue, name in TISSUES.items():
        #     print name, len(vhl_muts[vhl_muts["tissue"] == tissue]["sample"].unique())
        #print vhl_muts.to_csv(sep="\t", index=False),

        print results_df.to_csv(sep="\t", index=False),
