# Time-Domain-THz-Analysis-Code

This is a comprehensive collection of data analysis routines for analyzing time-domain terahertz spectroscopy (TDTS) data with [Igor Pro](https://www.wavemetrics.com/products/igorpro/igorpro.htm), a data analysis application provided by [WaveMetrics](https://www.wavemetrics.com/index.html).

IgorPro.  This code was originally developed for use in Dr. Peter Armitage's Complex Materials Spectroscopy Laboratory at The Johns Hopkins University.  However, the code is robust and should be easily adaptable by any experimentalist who wishes to analyze TDTS data.

Detailed instructions for using each piece of code can be found at the top of each .ipf file.  However, a brief summary of the code is as follows:

## **MAIN ANALYSIS PROCEDURE**:

The main analysis procedure is the "NJL_TDTSAnalysisCode.ipf" file.  This is the master code which is capable of loading raw data from a TDTS experiment, analyzing the code according to user specifications, and producing publishing quality plots.  In this code, the user specifies the type of sample being investigated - possibilities include single crystals, thin films, or single crystals mounted to a substrate - and also the type of data being taken - either as a function of temperature or magnetic field.  Although, the code should be easily adaptable to analyze data as a function of any parameter, it's just a matter of naming convention.  This procedure is also capable of analyzing data taken with the fast rotating polarizer technique developed in the Armitage lab. 

The user then loads the raw data via the analysis palette, which is accessible by holding CTRL-2, and then the code analyzes the data using the appropriate algorithms depending on the sample type and experiment settings.  The code is capable of computing the complex conductivity, the complex magnetic susceptibility, and Faraday rotation in both the linear and circular bases.  It can also perform frequency cuts at user specified frequencies for all the above quantities.

## **ONLY CALCULATING THE MAGNETIC SUSCEPTIBILITY OF A SAMPLE**:

I've also written a new procedure "NJL_MagneticSusceptibilityCode.ipf" which allows the user to quickly calculate the magnetic susceptibility of a sample.  This code is a greatly condensed and simplier version of the main analysis procedure described above.  Here, the user feeds the program only the complex raw transmission of a sample (with naming conventions described in the code) and then code calculates the complex index of refraction and complex magnetic susceptibility of the sample.  This code is very user friendly and should make calculating the magnetic susceptibilty of a sample quick and easy.

## **USING THIS CODE**:

I wrote this code in the earlier days of my Ph.D. to help streamline the data analysis process for both myself and my lab-mates.  My genuine hope is that it will be useful to the scientific community at large.  With that being said, I ask that you please give credit if using or adapting this code.  I request that you acknowledge the use of this code in any manuscript for which it was used as: "We thank N. J. Laurita for access to his data analysis routines."  Also, If you compute the magnetic susceptibility of your sample with these procedures, then please cite my PRL from which that analysis originated: https://link.aps.org/doi/10.1103/PhysRevLett.114.207201

If you found this code useful or have suggestions for improvements please let me know! I'd love to hear about the exciting uses that this code has contributed to, I can be reached at NLaurita1@Gmail.com.

Best of luck-
NJL
