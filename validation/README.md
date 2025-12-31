# Validation Folder

Scripts for validating toolkit predictions against known data.

## Validation Tests:

1. **validate_lipinski.py** - Test RO5 filter accuracy
2. **benchmark_fda_drugs.py** - Compare predictions with FDA-approved drugs

## Running Validation:
```bash
python validation/validate_lipinski.py
python validation/benchmark_fda_drugs.py
```

All validation scripts generate detailed reports and statistics.