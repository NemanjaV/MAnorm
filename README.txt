1. Pre-requisites:
     a. Bedtools installed: http://code.google.com/p/bedtools/                                                                                                                                                      
     b. Bioconductor packages installed: affy,R.basic and MASS                                                                                                                                                      
       HOWTO: input the following lines to install these 3 packages in R                                                                                                                                            
            #affy:                                                                                                                                                                                                  
                  source("http://bioconductor.org/biocLite.R")                                                                                                                                                      
                  biocLite(c("affy")                                                                                                                                                                                
            #R.basic and MASS:                                                                                                                                                                                      
                     source("http://www.braju.com/R/hbLite.R")                                                                                                                                                      
                     hbLite(c("R.basic","MASS"))                                                                                                                                                                    
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
2. Usage:
       I. Create a folder and place in the folder MAnorm.sh, MAnorm.r, and all 4 bed files (peak file and read file for the 2 samples under comparison) to be analyzed.

      II. run command:   ./MAnorm.sh   sample1_peakfile[BED]     sample2_peakfile[BED]     sample1_readfile[BED]    sample2_readfile[BED]     sample1_readshift_lentgh[INT]       sample2_readshift_length[INT]                        
                                                                                                                                                                                                                    
3. Example command:                                                                                                                                                                                                 
      ./MAnorm.sh sample1_peaks.bed sample2_peaks.bed sample1_read.bed sample2_read.bed 150  150                                                                                                                     
                                                                                                                                                                                                                    
      The first 4 parameters should be input files in bed format with no header lines
         The first 2 files have ONLY 3 columns: chromosome, start, end.
         The next 2 files should have 4 columns: chromosome, start, end, strand (+/-)

      The last 2 parameters are the number of bp to be shifted for each read. These two parameters are found from MACS peak file *_peaks.xls after "# d =".
                                                                                                                                                                                                                    
      MAnorm.r is called from MAnorm.sh, and there is no need to run it separately.  It checks file Rcommand.out for the output file from running the R script.                                                                       
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
4. Output:                                                                                                                                                                                                          
                                                                                                                                                                                                                    
4.1 Output files created:                                                                                                                                                                                           
                                                                                                                                                                                                                    
     MAnorm_result.xls: Includes all peak coordinates, # raw reads from sample1, #raw reads from sample2, M-value_rescaled, A-value_rescaled, -log10(p-value). 
     This file lists the peak coordinates, raw reads, M and A values, and MA-norm p-value for the set of common peaks in each sample and for the peaks unique to each sample.  
     The file below (MAnorm_result_commonPeak_merged.xls) provides the list of merged common peaks.                                                  
           
     MAnorm_result_commonPeak_merged.xls: This output file is similar to MAnorm_result.xls, except that the common peaks are merged. This file corresponds to a list of merged peaks for the two samples being compared.
                                                                                                       
     sample1_peaks.wig: wig file for sample1 peak list, with the width of the bar representing the M value                                                                                                   
     sample2_peaks.wig: wig file for sample2 peak list, with the width of the bar representing the M value                                                                                                   
                                                                                                                                                                                                                    
     MAplot_before_rescaling (common peaks).png : MAplot for common peaks before MAnorm                                                                                                                      
                   green line: robust regression                                                                                                                                                                    
                   red line: LOWESS regression                                                                                                                                                                       
                   blue line: line M = 0                                                                                                                                                                            
     MAplot_after_rescaling (all peaks).png : MAplot for all peaks after MAnorm                                                                                                                              
                   red line: LOWESS regression                                                                                                                                                                       
                   blue line: line M = 0                                                                                                                                                                            
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
4.2 Temporary Output files removed after MAnorm is done(If these files are needed, please modify MAnorm.sh):                                                                                                                                                                  
                peak1.bed: sample 1 peak list
                peak2.bed: sample 2 peak list
                read1.bed: sample 1 read used for mapping to peaks
                read2.bed: sample 2 read used for mapping to peaks
                peak1_dump.bed: sample 1 peaks not used
                peak2_dump.bed: sample 2 peaks not used
                read1_dump.bed: sample 1 read NOT used for mapping to peaks
                read2_dump.bed: sample 2 read NOT used for mapping to peaks                                                                                                                                              
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
5.  classfy_by_M_value.sh is a useful script  for peak classification based on M value, which could be uploaded directly to a genome browser.  Required to be in the same folder of MAnorm_result.xls generated by MAnorm.sh.
                                                                                                                                                                                                                    
     usage:                                                                                                                                                                                                         
               ./classfy_by_M_value.sh       absolute_M-value_cutoff_unbiased          absolute_M-value_minimum_for_unique_peaks               -log10p_cutoff_unique
                                                                                                                                                                                                                    
     example_command:                                                                                                                                                                                               
               ./classfy_by_M_value.sh 0.5  1  5

In the above example, sample 1 unique peaks (non-concordant peaks) are peaks with M-values greater than 1 and that have a log base 10(p-value) greater than 5.  Similarly, sample 2 unique peaks (non-concordant peaks) are peaks with M-values less than (-1) that have a log base 10(p-value) greater than 5.  Unbiased peaks (concordant peaks) are peaks with M-values between (-0.5) and (+0.5).                                                                                                                                                                       
                                                                                                                                                                                                                    
     output:                                                                                                                                                                                                        
               3 bed files:  2 unique peak lists and 1 unbiased peak list (Example of output files in “Output” folder.)                                                                                                                                            
                                                                                                                                                                                                                    
6. contact:                                                                                                                                                                                                         
       zhangyij@gmail.com                                                                                                                                                                                           
