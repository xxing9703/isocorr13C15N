# isocorr13C15N -- isotope correction for both 13C and 15N

isocorr13C15N a matlab gui version of the natural isotope correction tool that supports both 13C and 15N labeled experiments. It also works for single 13C or 15N labeled data inputs.  The algorithm is based on the following paper and associated R code:

Su X, Lu W and Rabinowitz J (2017). "Metabolite Spectral Accuracy on Orbitraps." Analytical Chemistry, 89(11), pp. 5940-5948. doi: 10.1021/acs.analchem.7b00396 (URL: http://doi.org/10.1021/acs.analchem.7b0039), PMID: 28471646, R package version 0.2.3 (2018).
R code:  https://github.com/XiaoyangSu/AccuCor

It is simplified to only consider high resolution mass spectrometer that can well resolve 13C and 15N, which is the default for Exactive, commonly used in the lab. Therefore, instrumental resolution is no longer an input parameter that users need to enter. C and N purity is 99% as the default, which can be changaed from the interface.

# Input file

 ".csv" file directly from Maven or el-Maven output is supported, which includes the following key columns:
 -metaGroupId
 -isotopeLabel
 -compound
 -formula
 Data block always start from column 15
 
 The isotope for correction is automatically determined by the contents of the 'isotopeLabel' column.
 The formula name string in the formula column needs to be correct in order to determin the C/N numbers.
 
 Example input files are included, where example1.csv is 13C only, and example2.csv contains both 13C and 15N.
 
 # Output file
 
 Output file name for "input.csv" will be "input_cor.csv", automatically created in the same folder.
 
 # Enviorment
  requires Matlab 2018b or later being installed
  
 Follow instructions on the GUI.


 # previous update:

7/17 use metaGroupID as identifier,instead of formula name, to group data blocks belonging to the same metabolite. Sometimes, two metabolites with the same names/formula but different rt are considered different entries.

8/16 readtable(~, 'variablename', true), allows column headers of the input file to be numbers or others. 

8/17 use 'try... catch' to parse formula string, provides warning messagebox with #row hint if formula of that row contains errors.
 

