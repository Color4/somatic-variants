#
# Makefile for setting up and running analyses.
#

# Set up a conda environment for analysis.
setup:
	conda create -n somatic-variants-env -y python
	source activate somatic-variants-env && \
		conda install -c bioconda -y pandas gemini && \
		pip install toyplot

# Compare variants in normal and serum.
serum_vs_normal: augmented_find_vars_results_all_minaf0.02_ad5_d20.txt
	source activate somatic-variants-env && \
	python ~/projects/somatic-variants/gemini/compare_tissue_variants.py augmented_find_vars_results_all_minaf0.02_ad5_d20.txt > serum_vs_normal_summary.txt

# Get variants meeting minimum filtering criteria that are in normal or serum.
augmented_find_vars_results_all_minaf0.02_ad5_d20.txt:
	source activate somatic-variants-env && \
	python ~/projects/somatic-variants/gemini/rcc.py augment_vars --results-file find_vars_results_all_minaf0.02_ad5_d20.txt --tissue N,S


# Remove analysis environment.
clean:
	conda remove -n somatic-variants-env -y --all
