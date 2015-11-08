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
from toyplot import Canvas, browser

# Columns for output.
DISPLAY_COLS = ["sample", "id", "plate", "tissue", "replicate"] + \
               ["gene", "chrom", "start", "ref", "alt", "num_het", "num_het_by_id", "impact", "impact_severity", "cosmic_ids", "rs_ids", \
                "sift_pred", "polyphen_pred", "vep_hgvsc", "vep_hgvsp"] + \
               ["alt_depth", "depth", "allele_freq"]


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
    NOTE: this funciton is for dataframes where each variant (not genotype) is a single row with genotypes information as columns.
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
    NOTE: this funciton is for dataframes where each variant (not genotype) is a single row with genotypes information as columns.
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


def get_variants_in_samples(gemini_db, samples, annotations, min_allele_freq, min_alt_depth, min_depth, max_aaf_all, somatic=False):
    """
    Returns a dataframe with variants from all samples.
    """

    if somatic:
        get_vars_fn = gem_ops.get_somatic_vars_in_sample
    else:
        get_vars_fn = gem_ops.get_vars_in_sample

    print >> sys.stderr, "Samples to process: ", ", ".join(samples)
    all_vars = []
    for sample in samples:
        start = time.time()
        sample_vars = get_vars_fn(gemini_db, annotations, sample, min_allele_freq=min_allele_freq,
                                  min_alt_depth=min_alt_depth, min_depth=min_depth, max_aaf_all=max_aaf_all)
        end = time.time()
        all_vars.append(sample_vars)
        print >> sys.stderr, sample, len(sample_vars), end-start

    # Combine all variants together and reset index so adding metadata is easy.
    all_vars_df = gem_ops.convert_cols(pd.concat(all_vars))
    all_vars_df.reset_index(inplace=True, drop=True)

    # Using sample column, create and populate columns for id, plate, tissue, and replicate.
    sample_attrs = pd.DataFrame(list(all_vars_df["sample"].apply(lambda s: split_id(s))), columns=["id", "plate", "tissue", "replicate"])
    for i, col in enumerate(["id", "plate", "tissue", "replicate"]):
        all_vars_df.insert(i+1, col, pd.Series())
        all_vars_df[col] = sample_attrs[col]

    return all_vars_df

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

    column_templates = {
        "sn": "j_Serum%i_Normal%i",
        "st": "j_Serum%i_Tumor%i",
        "tn": "j_Tumor%i_Normal%i",
        "stn": "j_Serum%i_Tumor%i_Normal%i",
        "sf": "j_Serum%i_FFPE",
        "nf": "j_Normal%i_FFPE",
        "tf": "j_Tumor%i_FFPE",
        "stnf": "j_Serum%i_Tumor%i_Normal%i_FFPE"
    }

    # Add columns for joint variants.
    # NOTE: should set max serum/tumor/normal by looking at data frame and grouping
    # by ID + tissue.
    num_serum_samples = 2
    num_tumor_samples = 2
    num_normal_samples = 2
    for s in range(num_serum_samples):
        for t in range(num_tumor_samples):
            for n in range(num_normal_samples):
                sn_col = column_templates["sn"] % (s, n)
                st_col = column_templates["st"] % (s, t)
                tn_col = column_templates["tn"] % (t, n)
                stn_col = column_templates["stn"] % (s, t, n)
                sf_col = column_templates["sf"] % (s)
                nf_col = column_templates["nf"] % (n)
                tf_col = column_templates["tf"] % (t)
                stnf_col = column_templates["stnf"] % (s, t, n)
                somatic_vars[sn_col] = False
                somatic_vars[st_col] = False
                somatic_vars[tn_col] = False
                somatic_vars[stn_col] = False
                somatic_vars[sf_col] = False
                somatic_vars[nf_col] = False
                somatic_vars[tf_col] = False
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
                    sn_col = column_templates["sn"] % (s, n)
                    somatic_vars.loc[var_group.index, sn_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, normal_sample])

                    # S and F
                    sf_col = column_templates["sf"] % (s)
                    somatic_vars.loc[var_group.index, sf_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample], any_has_gt_samples=ffpe_samples)

                    # N and F
                    nf_col = column_templates["nf"] % (n)
                    somatic_vars.loc[var_group.index, nf_col] = joint_genotypes(var_group, all_have_gt_samples=[normal_sample], any_has_gt_samples=ffpe_samples)

                    for t, tumor_sample in enumerate(tumor_samples):
                        # T and N
                        # T and F
                        # S and T
                        # S, T, and N
                        # S, T, N, and F
                        tn_col = column_templates["tn"] % (t, n)
                        tf_col = column_templates["tf"] % (t)
                        st_col = column_templates["st"] % (s, t)
                        stn_col = column_templates["stn"] % (s, t, n)
                        stnf_col = column_templates["stnf"] % (s, t, n)
                        somatic_vars.loc[var_group.index, st_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, tumor_sample])
                        somatic_vars.loc[var_group.index, tf_col] = joint_genotypes(var_group, all_have_gt_samples=[tumor_sample], any_has_gt_samples=ffpe_samples)
                        somatic_vars.loc[var_group.index, tn_col] = joint_genotypes(var_group, all_have_gt_samples=[tumor_sample, normal_sample])
                        somatic_vars.loc[var_group.index, stn_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, normal_sample, tumor_sample])
                        somatic_vars.loc[var_group.index, stnf_col] = joint_genotypes(var_group, all_have_gt_samples=[serum_sample, tumor_sample, normal_sample], any_has_gt_samples=ffpe_samples)

                if len(serum_samples) == 0:
                    for t, tumor_sample in enumerate(tumor_samples):
                        tn_col = column_templates["tn"] % (t, n)
                        somatic_vars.loc[var_group.index, tn_col] = joint_genotypes(var_group, all_have_gt_samples=[tumor_sample, normal_sample])

    return somatic_vars

def filter_and_augment_variants(vars_df, min_allele_freq, min_alt_depth, min_depth, max_num_het, tissue, add_joint_cols=True):
    """
    Filter and augment variants dataframe.
    """

    # Filter using allele freq, alt depth, and depth.
    vars_df = vars_df[(vars_df.allele_freq >= min_allele_freq) &
                      (vars_df.alt_depth >= min_alt_depth) &
                      (vars_df.depth >= min_depth)]

    # Update num_het, add num_het_by_id, germline
    vars_df = gem_ops.update_num_het(vars_df)
    vars_df = gem_ops.update_num_het_by_id(vars_df)
    vars_df = gem_ops.update_germline(vars_df)
    vars_df = gem_ops.update_mean_af(vars_df)

    # Filter by max_num_het and tissue.
    vars_df = vars_df[(vars_df.num_het <= max_num_het) & (vars_df.tissue).isin(tissue)]

    # Add joint columns.
    if add_joint_cols:
        vars_df = add_joint_cols(vars_df)

    return vars_df

def print_joint_variants(somatic_vars, stream, all_cols=False):
    """
    Print joint variants.
    """
    if all_cols:
        cols = somatic_vars.columns.values
    else:
        cols = DISPLAY_COLS

    for an_id, id_vars in somatic_vars.groupby("id"):
        print >> stream, "\n%s ****************" % an_id

        print >> stream, "Number of variants per sample:"
        print >> stream, id_vars.groupby("sample").size()

        joint_cols = [col for col in somatic_vars.columns.values if col.startswith("j_")]
        for col in joint_cols:
            col_tissues = map(lambda s: s[0], col.split("_")[1:])
            id_shared_vars = id_vars[(id_vars[col] == True) & id_vars["tissue"].isin(col_tissues)]

            print >> stream, "\n%s: %i shared variants ----------" % (" + ".join(col.split("_")[1:]), len(id_shared_vars["variant_id"].unique()))
            if len(id_shared_vars) > 0:
                #print >> stream, "\n%s: %i shared variants ----------" % (" + ".join(col.split("_")[1:]), len(id_shared_vars["variant_id"].unique()))
                print >> stream, id_shared_vars.sort("variant_id")[cols].to_csv(sep="\t", index=False),

        # Print variant in serum but not in normal.
        print >> stream, "\nVariants in Serum but not Normal ----------"
        for variant, genotypes in id_vars.groupby("variant_id"):
            if len(genotypes[genotypes.tissue == "S"]) > 0 and len(genotypes[genotypes.tissue == "N"]) == 0:
                print >> stream, genotypes[genotypes.tissue == "S"][cols].to_csv(sep="\t", index=False, header=False),

        print >> stream, "\nVariants in Tumor but not Normal ----------"
        for variant, genotypes in id_vars.groupby("variant_id"):
            if len(genotypes[genotypes.tissue == "T"]) > 0 and len(genotypes[genotypes.tissue == "N"]) == 0:
                print >> stream, genotypes[genotypes.tissue == "T"][cols].to_csv(sep="\t", index=False, header=False),


if __name__ == "__main__":
    # Argument setup and parsing.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("operation", help="Operation to perform")
    parser.add_argument("--gemini-db", help="Gemini database to use")
    parser.add_argument("--min-allele-freq", type=float, default=0.02, help="Allele frequency to use")
    parser.add_argument("--max-aaf-all", type=float, default=0.01, help="Maximum alternate allele frequency")
    parser.add_argument("--min-depth", type=int, default=20, help="Minimum depth to use")
    parser.add_argument("--min-alt-depth", type=int, default=5, help="Minimum alternate depth to use")
    parser.add_argument("--max-num-het", type=int, default=sys.maxint, help="Maximum number of observations")
    parser.add_argument("--tissue", default="S,T,N,F", help="Filter for only these tissue")
    parser.add_argument("--sample-pattern", default=".*", help="Regular expression to use to select samples")
    parser.add_argument("--annotations", default="TCGA_RCC", help="Annotations to use as hotspots")
    parser.add_argument("--results-file", default="", help="File of raw somatic results")
    parser.add_argument("--add-joint", action="store_true", help="File of raw somatic results")

    args = parser.parse_args()

    gemini_db = args.gemini_db
    operation = args.operation
    min_allele_freq = args.min_allele_freq
    min_depth = args.min_depth
    min_alt_depth = args.min_alt_depth
    max_num_het = args.max_num_het
    max_aaf_all = args.max_aaf_all
    tissue = args.tissue.split(",")
    sample_pattern = args.sample_pattern
    annotations = args.annotations.split(",")
    results_file = args.results_file
    add_joint = args.add_joint

    if operation == "find_all":
        # Get all variants in a set of samples.

        # Get samples to process.
        samples = [sample for sample in gem_ops.get_samples(gemini_db) if re.search(sample_pattern, sample) > 0]

        # Get sample variants.
        all_vars_df = get_variants_in_samples(gemini_db, samples, annotations, min_allele_freq, min_alt_depth, min_depth, max_aaf_all, somatic=False)

        # Write results to file.
        if sample_pattern == ".*":
            sample_pattern = "all"
        out_filename = "find_vars_results_%s_minaf%.2f_ad%i_d%i.txt" % (sample_pattern, min_allele_freq, min_alt_depth, min_depth)
        out_file = open(out_filename, "w")
        out_file.write( all_vars_df.to_csv(sep="\t", index=False, float_format='%.3f') )
        out_file.close()

        print "Wrote results to file %s" % out_filename

    elif operation == "augment_vars":
        # Augment variants with updated and joint information.

        # Read results into dataframe.
        results_df = gem_ops.convert_cols( pd.read_csv(results_file, sep="\t") )
        results_df = filter_and_augment_variants(results_df, min_allele_freq, min_alt_depth, min_depth, max_num_het, tissue, add_joint_cols)

        # Print augmented results.
        augmented_out_file = open("augmented_" + results_file, "w")
        augmented_out_file.write(results_df.to_csv(sep="\t", index=False))
        augmented_out_file.close()

        print "Wrote augmented variants to file %s" % augmented_out_file

        # Print joint variants.
        if add_joint:
            joint_out_file = open("joint_" + results_file, "w")
            print_joint_variants(results_df, joint_out_file, all_cols=True)
            joint_out_file.close()

        print "Wrote joint variants to file %s" % joint_out_file

    elif operation == "sweeps":
        # Perform sweeps using AF and DP/AD.

        # Read results into dataframe.
        results_df = gem_ops.convert_cols( pd.read_csv(results_file, sep="\t") )

        # Remove artifacts and germline variants.
        results_df = results_df[(results_df.HP <= 4) & (results_df.num_het <= 10) & (results_df.germline == False)]

        # Sweep across allele frequency and depths and compute number of variants and number of samples showing a variant
        # for each AF and DP.
        genes = ["VHL", "PBRM1", "MUC4", "SETD2", "BAP1", "KDM5C"]
        allele_freqs = [af/100.0 for af in range(2, 20)]
        depths = range(10, 200, 30)
        num_samples = len(results_df["sample"].unique())

        for gene in genes:
            canvas = Canvas(width=500, height=500)
            axes = canvas.axes(label=gene, ymin=0)

            gene_df = results_df[results_df.gene == gene]
            for d in depths:
                x = []
                y = []
                s = []
                for af in allele_freqs:
                    kept_vars = gene_df[(gene_df.allele_freq >= af) &
                                        (gene_df.alt_depth >= d) &
                                        (gene_df.impact_severity.isin(["MED", "HIGH"]))]
                    print "%s\t%0.2f\t%i\t%i\t%i\t%0.2f" % (gene, af, d, len(kept_vars), len(kept_vars["sample"].unique()), len(kept_vars["sample"].unique())/float(num_samples))
                    x.append(af)
                    y.append(d)
                    s.append(len(kept_vars["sample"].unique())/float(num_samples))
                mark = axes.plot(x, s)
            browser.show(canvas)

    elif operation == "augmented_to_somatic":
        # Go from augmentd variants to somatic.

        # TODO: not finished.
        pass

        # # Read results into dataframe.
        # results_df = gem_ops.convert_cols( pd.read_csv(results_file, sep="\t") )
        #
        # # Apply filtering criteria:
        # results_df = results_df[(results_df.HP < 5) & (results_df.rmsk == "None") & results_df.allele_freq >= 0.02 ]


    # elif operation == "find_somatic":
    #     # Get somatic variants at the individual sample level.
    #
    #     # Get samples to process.
    #     samples = [sample for sample in gem_ops.get_samples(gemini_db) if re.search(sample_pattern, sample) > 0]
    #
    #
    #     all_somatic_df = get_variants_in_samples(gemini_db, samples, annotations, min_allele_freq, min_alt_depth, min_depth, somatic=True)
    #
    #     # Write results to file.
    #     if sample_pattern == ".*":
    #         sample_pattern = "all"
    #     out_filename = "find_somatic_results_%s_minaf%.2f_ad%i_d%i.txt" % (sample_pattern, min_allele_freq, min_alt_depth, min_depth)
    #     out_file = open(out_filename, "w")
    #     out_file.write( all_somatic_df.to_csv(sep="\t", index=False, float_format='%.3f') )
    #     out_file.close()
