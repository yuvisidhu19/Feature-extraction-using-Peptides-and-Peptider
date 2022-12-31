# Feature extraction using Peptides and Peptider libraries
This R script can be used to extract multiple features from the peptide sequence provided in input.csv file using peptides and peptider libraries. It generates an output.csv file with the new features and a log file.

This program is capable to handle exception/error (if any) and write it to log file:
- Correct number of parameters (inputFileName).
- Show the appropriate message for wrong inputs.
- Handling of “File not Found” exception
- Many sequences are missing in the input file, handle them (ignore them).
- If any issue with the input record, it must be write to a log file

This program can only be run from command line using the following syntax: rscript <program.r> <InputDataFile>

E.g., ```rscript feature_extraction.r input.csv```
