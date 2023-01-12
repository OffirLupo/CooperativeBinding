# CooperativeBinding
Code for analysis and visualization of data presented in the manuscript "The architecture of binding cooperativity between densely bound transcription factors".

This code runs on out files (genome tracks of aligned reads) found in GEO accession  ####.

First step for all figures is to load a matlab struct will all out files loaded and normalized - using the function chec_struct_paired.m to load paired-end out files
and normalize the data for read count. second step is averaging repeats using the function mean_of_repeats2.m.
This processed matlab struct will be available at dyrad soon.
This mat struct needs to be loacted in a subfolder (e.g, Data_structs/checData.m)

Several functions in this repository, e.g., cbrewer.m and  joyPlot.m were taken from Matlab File Exchange. 
We thank the authors for sharing these files with the community, and ask them to contact if there are any issues with including those files in our repository. 

If experiencing any issues, or have any other comments/question, please contact (offirlupo1@gmail.com) 

Last update: 12/01/2023
