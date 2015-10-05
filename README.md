Processing variants to find somatic variants without a matched normal.
================

Creating a complete GEMINI database
----------------

1. Modify settings/paths in `gemini/default_settings`
2. Make sure requirements in `gemini/create_gemini_db.sh`
3. Run `gemini/create_gemini_db.sh`

Querying an RCC GEMINI database
----------------

NOTE: for script to run, your PYTHONPATH must include GEMINI; see `gemini/default_settings.sh` for an example.

1. Modify settings/paths in `gemini/default_settings`
2. Run `PYTHONPATH=/path/to/gemini python /this/directory/gemini/rcc.py` to see options for filtering variants
3. Example usage: `PYTHONPATH=/path/to/gemini python /this/directory/gemini/rcc.py --min-allele-freq 0.05 --min-alt-depth 5 --min-depth 20 --sample-pattern tumor my_variants.db`
