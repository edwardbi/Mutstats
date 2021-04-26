The Mutstats method takes the following arguments:
1) sample_id: the name or id of the sample
2) input_file_dir: the name of the input file of the sample
3) R_code_sir_dir: the code directory of mutstat.R
4) output_dir: the folder directory to contain the output

When running Mutstats on a sample, type in the following commend:
R CMD BATCH --no-save --no-restore '--args sample_id input_file_dir R_code_sir_dir output_dir' R_code_sir_dir/mutstat.R ./sample.log