Processing variants to find somatic variants without a matched normal.
================

Creating a complete GEMINI database
----------------

1. Modify settings/paths in `gemini/default_settings`
2. Make sure requirements in `gemini/create_gemini_db.sh`
3. Run `gemini/create_gemini_db.sh`

Querying an RCC GEMINI database for somatic variants.
----------------

NOTE: for script to run, your PYTHONPATH must include GEMINI; see `gemini/default_settings.sh` for an example.

1. Modify settings/paths in `gemini/default_settings`
2. Run `PYTHONPATH=/path/to/gemini python /this/directory/gemini/rcc.py` to see options for filtering variants
3. Example usage: `PYTHONPATH=/path/to/gemini python /this/directory/gemini/rcc.py find_somatic --gemini-db my_variants.db --min-allele-freq 0.05 --min-alt-depth 5 --min-depth 20 --sample-pattern tumor > somatic.txt`
4. Augment results with joint columns and num_het_by_id, fix num_het: `PYTHONPATH=/path/to/gemini python /this/directory/gemini/rcc.py augment_somatic --results-file somatic.txt > augmented_somatic.txt`
