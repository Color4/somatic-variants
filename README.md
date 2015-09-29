Processing variants to find somatic variants without a matched normal.
================

Creating a complete GEMINI database
----------------

1. Modify settings/paths in `gemini/default_settings`
2. Make sure requirements in `gemini/create_gemini_db.sh`
3. Run `gemini/create_gemini_db.sh`

Querying an RCC GEMINI database
----------------

1. Modify settings/paths in `gemini/default_settings`
2. Run `python /this/directory/gemini/rcc.py` to see options for filtering variants
3. Example usage: `python /this/directory/gemini/rcc.py --min-allele-freq 0.01 --min-alt-depth 1 --min-depth 1 my_variants.db > output.txt`
