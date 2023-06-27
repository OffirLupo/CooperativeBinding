# CooperativeBinding
Code for analysis and visualization of data presented in the manuscript "The architecture of binding cooperativity between densely bound transcription factors".

***to run with .out files:*** 

out files are genome tracks of aligned reads found in GEO accession GSE222268.
First step for all figures is to load a matlab struct will all out files loaded and normalized - using the function chec_struct_paired.m to load paired-end .out files
and normalize the data for read count.
second step is averaging repeats using the function mean_of_repeats2.m.

**These norm files (average of repeats) are also found in folders normFiles and WT_normFiles:**

genomic tracks (normalized for total read count) are found in this repository and can be loaded directly using the function loadMatFiles.m.
save this combined mat struct in a subfolder (e.g, Data_structs/checData.mat)

Several functions in this repository, e.g., cbrewer.m and  joyPlot.m were taken from Matlab File Exchange. 
We thank the authors for sharing these files with the community, and ask them to contact if there are any issues with including those files in our repository. 

If experiencing any issues, or have any other comments/question, please contact (offirlupo1@gmail.com) 

Last update: 15/05/2023
