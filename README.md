# Time-Domain-THz-Analysis-Code

This is a comprehensive collection of data analysis routines for analyzing time-domain terahertz spectroscopy (TDTS) data with IgorPro.  This code was originally developed for use in Dr. Peter Armitage's Complex Materials Spectroscopy Laboratory at The Johns Hopkins University.  However, the code is robust and should be easily adaptable by any experimentalist who wishes to analyze TDTS data.

Detailed instructions for using each piece of code can be found at the top of the .ipf file.  However, a brief summary of the code is as follows:

The main code is the "NJL_TDTSAnalysisCode.ipf" file.  This is the master code which is capable of loading raw data from a TDTS experiment, analyzing the code according to use specifications, and producing publishing quality plots.  In this code, the user specifies the type of sample of which spectroscopy is performed, possibilities include single crystals, thin films, or single crystals mounted to a substrate, and also the type of data being taken, either temperature dependent or magnetic field dependent.  However, the code should be easily adaptable to analyze data as a function of any parameter.  The code is also capable of analyzing data taken with the fast rotating polarizer technique developed in the Armitage lab. 

The user then loads the raw data via the analysis palette and then the code analyzes the data using the appropriate algorithms depending on the sample type and experiment settings.  The code is capable of computing the complex conductivity, the complex magnetic susceptibility, and Faraday rotation in both the linear and circular.  The code can also create folders with frequency cuts at the users specified frequencies of all the above quantities.

USING THIS CODE:
I wrote this code in the earlier days of my Ph.D. to help streamline the data analysis process for both myself and my lab mates.  My genuine hope is that by posting it here it will be useful to the TDTS community at large.  With that being said, I ask that you please give credit if using or adapting this code.  I request a simple acknowledgement in any publication from which this code has been used: "We thank N. J. Laurita for access to his data analysis routines."  Also, If you compute the magnetic susceptibility of your sample with this code then please cite my PRL from which that analysis originated: DOI: https://doi.org/10.1103/PhysRevLett.115.019901

If you found this code useful or have suggestions for improvements please let me know! I'd love to hear about the exciting uses that this code has contributed to, I can be reached at NLaurita1@Gmail.com.

Best of luck,
NJL
