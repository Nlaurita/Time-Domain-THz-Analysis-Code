# Time-Domain-THz-Analysis-Code

This is a comprehensive collection of data analysis routines for analyzing time-domain terahertz spectroscopy (TDTS) data with [Igor Pro](https://www.wavemetrics.com/products/igorpro/igorpro.htm), a data analysis application provided by [WaveMetrics](https://www.wavemetrics.com/index.html).

## **Running The Procedures**
In order to run these procedures I recomment saving the .ipf files in the "User Files Directory" of Igor so that the procedures can be accessed by any Igor experiment.  Each .ipf file has a control panel GUI which makes running the analysis routines easy.  The GUI can be accessed either through the macros menu or simple by pressing Ctrl-2.  Detailed instructions for using each piece of code can be found at the top of each .ipf file.  

## **Description of Routines**

Detailed instructions for using each piece of code can be found at the top of each .ipf file.  However, a brief summary of the code is as follows:

### **Prerequisite Routines For Main Analysis Procedure**
The main analysis procedure "NJL_TDTSAnalysisCode.ipf" was written to interface with a data processing procedure written by [laserstonewall](https://github.com/laserstonewall?tab=repositories).  I recommend saving [this](https://github.com/laserstonewall/THzTDSProcedures/blob/master/THz_Procedures_02-12-2014.ipf) procedure in the "User Files Directory" along with the routines provided here for best results and improved capabiliities.

### **Main Analysis Procedure**:

The main analysis procedure, "NJL_TDTSAnalysisCode.ipf" file, is a data analysis routine which is capable of loading raw data from a TDTS experiment, performing analysis according to user specifications, and producing publishable quality plots.  In this routine, the user specifies the type of sample being investigated - possibilities include single crystals, thin films, or single crystals mounted to a substrate - and also if the experiment was performed as a function of either temperature or magnetic field.  Although, the code should be easily adaptable to analyze data as a function of any parameter, it's just a matter of naming conventions.  This procedure is also capable of analyzing data taken with the [fast rotating polarizer technique](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-20-11-12303) developed in the Armitage lab. 

The routine is applicable to dielectric conductors, magnetic insulators, and superconducting samples and is capable of computing the complex conductivity, the complex magnetic susceptibility, and Faraday rotation in both the linear and circular bases.  It can also perform frequency cuts at user specified frequencies for all the above quantities.

### **Calculating The Optical Complex Magnetic Susceptibility Of A Sample**:

Also included is a procedure "NJL_MagneticSusceptibilityCode.ipf" which allows the user to quickly calculate the complex magnetic susceptibility of a sample from a TDTS measurement.  This procedure is a greatly condensed and simplified version of the main analysis procedure described above.  Here, the user feeds the program only the complex raw transmission of a sample (with naming conventions described in the code) and then the procedure calculates the complex index of refraction and complex magnetic susceptibility of the sample.  For more information regarding the method for calculating the complex magnetic susceptibility, please see chatper 2 of my [PH.D. thesis](https://drive.google.com/file/d/0B7K_8wnuzg_QRE9lQkczZlpEUDQ/view).

## **Code Use**:

I wrote this code in the earlier days of my Ph.D. to help streamline the data analysis process for both myself and my lab-mates.  My genuine hope is that it will be useful to the scientific community at large.  With that being said, I ask that you please give credit if using or adapting this code.  I request that you acknowledge the use of this code in any manuscript for which it was used as: "We thank N. J. Laurita for access to his data analysis routines."  Also, If you compute the magnetic susceptibility of your sample with these procedures, then please cite my PRL from which that analysis originated: https://link.aps.org/doi/10.1103/PhysRevLett.114.207201

## **Credits**
These procedures were developed in the [group](https://sites.google.com/site/nparmitagegroup/) of Professor [N. Peter Armitage](http://physics-astronomy.jhu.edu/directory/n-peter-armitage/) in the [Johns Hopkins University Physics Department](http://physics-astronomy.jhu.edu/).

## **Contact**
If you would like help using this code, found this code useful, or have suggestions for improvements please feel free to contact me at Laurita.Nicholas@gmail.com.

