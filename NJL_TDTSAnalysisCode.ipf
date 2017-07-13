#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>

//************************************************************************************ PLEASE GIVE PROPER ACKNOWLEDGEMENT WHEN USING THIS CODE ************************************************************************************//
//
//  This code is public and it is my sincere hope that it will be useful to others in the field, in fact that's partially why I wrote it!  However, if this code does aid in the analysis of your data then I kindly ask that proper 
//  acknowledgement be given in return. 
//
//  I would greatly appreciate if you would do the following:
//  1) Acknowledge me in the ``Acknowledgements" section of your paper as "We thank N. J. Laurita for access to his data analysis routines."
//  2) If you compute the magnetic susceptibility of your sample with this code then please cite my PRL from which that analysis originated: https://link.aps.org/doi/10.1103/PhysRevLett.114.207201
//  3) Let me know! I can be reached at Laurita.Nicholas@gmail.com.  I'm interested to hear how and why the code was useful and (most of all) how it can be improved.
//
//*******************************************************************************************************************************************************************************************************************************************************************//



//*******************************************************************************************************************************************************************************************************************************************************************//
//
//  This is a procedure, written by Dr. Nicholas J. Laurita on Sept 22nd - 25th 2015, which analyzes both temperature and field dependent time-domain terahertz spectrsocopy data.  This code was originally designed for use in the 
//  Complex Materials Spectroscopy Laboratory of Dr. N. Peter Armitage at The Johns Hopkins University.  Much of the code is specific to the experiments at JHU.  However, aside from the functions which load the data below,
//  much of these routines could be easily adapted for analysis by any group.
//
//  The code is able to analyze data taken on both single crystal and thin film samples which can be either electric conductors, magnetic insulators, or superconductors.  It is also able to analyze polarization rotation data such as 
//  Faraday rotation in both the linear or circular reference frame.
//
//  Much more information on the analysis routines, particularly on computing the magnetic susceptibility in a time-domain THz measurement can be found in my Ph.D. thesis which can be found here:
//  https://drive.google.com/file/d/0B7K_8wnuzg_QRE9lQkczZlpEUDQ/view
//
//*******************************************************************************************************************************************************************************************************************************************************************//



//************************************************************************************************** INSTRUCTIONS FOR RUNNING THE CODE *************************************************************************************************************//
//
// #1	Select the type of experiment you have performed.  The code can analyze data taken as a function of temperature or magnetic field, although the choice only changes the type of naming convention the program looks 
//		for and could easily be adapted to analyze data as a function of any parameter.  You also have the option of including the rotating polarizer, i.e. the "rotator", in the analysis.  When this option is chosen the code
//		will analyze both the in phase and out of phase channels so that quantities such as Faraday rotation or data in the circular basis can be computed.  So say you performed an experiment as a function of magnetic
//		field with the rotator, you would then select both the "Magnetic Field" and "Use Rotator" check boxes on the panel.
//
// #2	Select the type of sample you are measuring.  Choices include single crystals with an aperture reference, thin films with a substrate reference, or single crystals with a substrate reference.  The choice here determines
//		which alogorithms are used to analyze the data.
//
// #3	Enter your sample name and reference name in the text box on the panel.  The conventions, which are slightly different for temperature dependent and magnetic field dependent data, are listed on the front panel.
//
// #4	The appropriate text waves must be created before hitting the "Get Folder Info" button.  These textwaves tell the code which folders to create for the analysis and where to load the data in the experiment. 
//
//				For temperature dependent data:
//					Create a text wave called "Temps" which contains the temperatures at which data was taken.  For instance elements in this wave might be "4K", "10K", etc.  Additionally create two numeric waves,
//					one called "Num_Ref_Scans" and another called "Num_Samp_Scans".  Populate these waves with the number of reference and samples scans taken at each temperature, the program
//					will use these values to average the scans together to make a single wave for higher signal to noise.  For instance, assume you measured a sample at 15K with 2 reference scans and 
//					3 sample scans.  Then there should be a row in the waves "Temps", "Num_Ref_Scans" and "Num_Samp_Scans" which reads "15K", "2", "3" across.
//
//				For magnetic field dependent data:
//					Create a text wave called "Temps" which contains the temperatures at which taken was taken.  For instance elements might be "4K", "10K", etc.  Additionally you must create a text wave
//					for each temperature that was measured with the name "Fields_Temps".  For instance if data was taken at 15K then there will be an additoinal text wave with the name "Fields_15K".
//					Populate this text wave with the magnetic fields which were measured at 15K.  The format requires two digits after the decimal point, so elements could be"1.50" or "15.00" etc.
//
// #5	Hit the "Get Folder Info" button to bring up a menu from which you should select any scan in the folder you wish to load data from.  In order for the code to work properly all data must be saved in the same folder.
//		After you select any scan the program will create the appropriate folders for your experiment.  The main folder will haev the name "SampleName_..._Analysis" depending on which type of data is being loaded.
//		Drag all the waves you created before hitting the "Get Folder Info" button into this folder.  This should include "Temps", "Num_Ref_Scans", and "Num_Samp_Scans" for temperature dependent data and "Temps"
//		and all the "Fields_Temps" waves for field dependent data.
// 
// #6	It is now time to load you THz data.  This program is meant to interface with an additional Igor program written by Dr. Chris M. Morris that manipulates and loads TDTS raw data.  Please see instructions for using
//		that palette contained within that code.  Once the settings on that palette are to your liking, select the appropriate "Load" button on the TDTS anlaysis palette.  This is fairly self explanatory, the "Load Reference
//		Scans In Phase" button will load all the reference data, with the settings selected on Chris's palette, into the "In_Phase" "time traces" folder and compute their FFTs.  The "Load Sample Scans In Phase" buttons 
//		does the same thing with the in phase sample data.  The "Out Of Phase" buttons are only needed if the "Use Rotator" check box is selected and are used when calculating Faraday rotation or analyzing data in the 
//		circular basis.
//
// #7	With all the raw data from the experiment now loaded, analysis is very fast and simple.  The "Calculate Transmissions" button takes the complex ratio of each sample FFT to the corresponding reference FFT to
//		calculate the complex transmission.  Folders are created for the transmissions and their magnitudes.  If temperature dependent data is being analyzed then each temperature will have a single transmission.  If magnetic
//		field dependent data is being anlayzed then every field at every temperature should have a transmission.
//
// #8 	The "Switch To Circular Basis" button can only be used if the "Use Rotator" checkbox was selected at the beginning of the analysis and if the "Out Of Phase" data was loaded.  This button takes the complex transmission
//		of the in phase data, the waves which begin "Txx_..." and the out of phase data "Tyx_...", and converts the data into right and left hand circular channels by T_L = Txx + i*Tyx and T_R = Txx - i*Tyx.
//
// #9	Calculating additional quantities, such as complex conductivity or magnetic susceptibility, requires additional information of the sample or reference.  For single crystal samples one must only imput the sample thickness in 
//		millimeters in the appropriate text box.  Thin films samples require both the substrates real part of the index of refraction, e.g. "3" or "2.5" etc., and the "delta L" value in mircons.  The delta L is the difference in thickness
//		between the reference substrate and the substrate the thin film was grown on.  Typical values will be between +/- 15 microns or so.  Single crystals mounted to a substrate require both the sample thickness and the 
//		substrate index of refraction.
//
// #10	With the addtional information included, the complex conducitivity is easily calculated by hitting the "Calculate Conductivity / Conductance" button.  Conductivity is calculatd for single crystal samples while conductance
//		is computed for thin films.  There is an additional checkbox labeled "Use Calculated Substrate Index Of Refraction" which is useful if the substrate's index of refraction has strong frequency dependence in the THz range.
//		When this checkbox is not selected then the code analyzes the data assuming the substrate's index of reftaction is a real constant given by the value inputted in the "Sub Index Of Refraction" textbox above.  A user should
//		select this textbox when they have performed a separate TDTS experiment on the substrate alone and found that the frequency dependence of the index of refraction is too large to ignore OR there is a significant
//		imaginary part of the index of refraction.  The user then inputs the name of the substrate and the code looks for the appropriate folders with the substrates name, grabs the complex index of refraction of the substrate, 
//		and uses it for calculating the complex conductance of the thin film.
//
// #11 	For magnetic insulators one is more interested in the magnetic susceptibility than the complex conductivity.  The "Calculate Magnetic Susceptibility" button calculates the magnetic susceptibility of the sample
//		as a function of temperature or magnetic field depending on the type of experiment.  However, in order to extract the magnetic susceptibility one must select a reference temperature and enter it into the 
//		"Reference Temp:" textbox.  An ideal reference temperature is a temperature at which the dielectric contributions to the index of refraction will not appreciably change below AND magnetic correlations have not yet begun.
//		I would suggest a temperature above a magnetic transition or above the appearance of a magnetic excitation but below roughly 100K.  For instance, you may chose 80K as a reference temperature and then should inpurt "80K"
//		into the text box.  For more information regarding why a reference temperature is needed please see Ch. 2 of my Ph.D. thesis.
//
// #12	The "Calculate Faraday Rotation" buttons calculates the complex Faraday rotation.  The "Use Rotator" checkbox must be selected at the beginning of the experiment to calculate the Faraday rotation and one must 
//		load the out of phase data.  If the "Switch To Circular Basis" button was selected earlier in the anlaysis then the Faraday rotation will computed in both the linear and circular bases.
//
// #13 	The "Generate Constant Field Folders" button is only applicable to experiments as a function of magnetic field.  This buttons takes the field dependent data, which is sorted by temperature by default in the program, 
//		and instead arranges it into folders labeled by field.  For instance, instead of having data organized  by all the fields measure at 10K, this button creates a folder at every specific field measured with all temperatures
//		at which that field was measured.
//
// #14	The "Frequency Cuts" section of the palette interpolates the data at the specified frequencies to generate folders of cuts at those frequencies of all relevant data.  The user inputs the starting, ending, and step size
//		frequencies, all in THz, and then the program will generate folders of frequency cuts at those frequencies for the quanties which have been previously calculated.  For instance, conductivity or magnetic susceptibility
//		etc.
//
// #15	The final button on the palette, the "Make Graph" button, makes a publishing ready plot of all the data located in the "active fiolder" in the data brower.  Its a great way to plot a lot of data very quickly.  The "Offset Amount"
//		option creates an offset between waves in the plot to better display the data.
//
//*******************************************************************************************************************************************************************************************************************************************************************//


//*******************************************************************************************************************************************************************************************************************************************************************//
// This function defines all our global variables and stores them in the folder "NickLTHzExtras," and well has generates and pulls up the analysis palette when the user presses Ctrl-2
//*******************************************************************************************************************************************************************************************************************************************************************//

//----You can use Control-2 to pull up the button which runs the code below----//
Menu "Macros"
	"Display Nick L THz Analysis Palette/2", NickTHz()
End

Function NickTHz()

	//Save the current data folder, b/c we will be moving to the new data folder to hold
	//the global variables, and we will need to go back to that
	String dfSave = GetDataFolder(1)
	
	//Create data folder in Packages, isolates it from the user so don't get a bunch of variables clogging stuff up
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S root:Packages:NickLTHzExtras

	//Brings window to front if open already, or opens if not
	DoWindow/HIDE=? $("NickL_THzAnalysis_Palette")
	if (V_flag != 0)
		DoWindow/F NickL_THzAnalysis_Palette;
	else
		Execute/Q "NickL_THzAnalysis_Palette()"
	endif
	
	//---Lets define some global variables---//

	//---Folder where data will be loaded from.  It gets generated by the "Got Folder Info" button in the analysis palette---//
	String/G base_folder
	
	//---Sample_Name stores the name to look for in the sample scans.  For instance in "SmB6_5K_1", Sample_Name = SmB6---//
	String/G Sample_Name = strVarOrDefault("root:Packages:NickLTHzExtras:Sample_Name", "nan") 
	
	//---Ref_Name stores the name to look for in the reference scans.  For instance in "Ap_5K_1", Ref_Name = Ap---//
	String/G Ref_Name = strVarOrDefault("root:Packages:NickLTHzExtras:Ref_Name", "nan") 
		
	//---Sub_Name stores the name of a substrate which has already been analyzed via the code---// 
	//---This is used to supply the complex index of refract of the substrate for calculating the complex conductance in the thin film case---//
	String/G Sub_Name = strVarOrDefault("root:Packages:NickLTHzExtras:Sub_Name", "nan") 
	
	//---dsc is the single crystal sample thickness, thickness has to be in mm's.---//
	variable/G dsc = NumVarOrDefault("root:Packages:NickLTHzExtras:dsc", 0) 
	
	//---n is the index of refraction of the substrate.---//
	variable/G ns = NumVarOrDefault("root:Packages:NickLTHzExtras:ns", 0) 
	
	//---DeltaL is the thickness difference between reference and sample substratein microns.---//
	variable/G DeltaL= NumVarOrDefault("root:Packages:NickLTHzExtras:DeltaL", 0) 
	
	//---Freq_Low_Limit is the starting frequency for the frequency cuts---//
	variable/G Freq_Low_Limit = NumVarOrDefault("root:Packages:NickLTHzExtras:Freq_Low_Limit", 0) 
	
	//---Cut resolution is the step size of the frequency cuts---//
	variable/G Cut_Resolution = NumVarOrDefault("root:Packages:NickLTHzExtras:Cut_Resolution", 0) 
	
	//---Freq_High_Limit is the ending frequency for the frequency cuts---//
	variable/G Freq_High_Limit = NumVarOrDefault("root:Packages:NickLTHzExtras:Freq_High_Limit", 0) 
			
	//---This string stores the measurment type so the conductivity can be calculated correctly---//
	String/G SampType = StrVarOrDefault("root:Packages:NickLTHzExtras:SampType", "nan") 
	
	//---This Variable Will Store the numeric value of the popup menu for sample type---//
	variable/G PopVal = NumVarOrDefault("root:Packages:NickLTHzExtras:PopVal", 0) 
	
	//---OSA is the offset amount entered into the analysis palette---//
	variable/G OSA = NumVarOrDefault("root:Packages:NickLTHzExtras:OSA", 0) 
	
	//---This string stores the offset type so we know which graph to apply the offset to---//
	String/G OsetType = StrVarOrDefault("root:Packages:NickLTHzExtras:OsetType", "nan") 
	
	//---RefTemp is the reference temperature used in the susceptibility calculation---//
	String/G RefTemp = StrVarOrDefault("root:Packages:NickLTHzExtras:RefTemp", "nan") 
	
	//---Go back to the original data folder---//
	SetDataFolder dfSave
	
	//---Let's open the "Make Traces Different" Panel so it doesn't throw an error when we try to use it later---//
	CreateKBColorizePanel()
	//---Hide the window for it though, we don't need to see it to use it---//
	DoWindow/hide=1 KBColorizePanel
	
End

//---Generates the THz Analysis Palette---//
Window NickL_THzAnalysis_Palette() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1084,160,1836,800)
	ModifyPanel frameStyle=3
	SetDrawLayer UserBack
	SetDrawEnv fsize= 16,fstyle= 1
	DrawText 289,22,"TDTS Analysis Palette"
	DrawText 7,191,"Sample Name Syntax: \"SampleName_50K_#"
	DrawText 7,209,"Reference Name Syntax: \"RefName_50K_#"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 382,65,"SIngle Crystals:"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 382,102,"Thin Films:"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 381,141,"Single Crystal On Substrate:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 377,188,750,188
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 509,449,"Frequency Cuts:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,371,375,371
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,595,750,595
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 224,621,"Please Read Instruction In Code Before Use"
	SetDrawEnv fsize= 14
	DrawText 214,637,"You can use CRTL-2 to pull up the analysis palette"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 375,523,750,523
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 116,105,"Select Sample Type:"
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 60,40,"Data Was Taken As A Function Of:"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 7,177,"Temperature Dependent Data:"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 7,233,"Magnetic Field Dependent Data:"
	DrawText 7,248,"Sample Name Syntax: \"50K_30.000_kG_SampleName\""
	DrawText 7,266,"Reference Name Syntax: \"50K_30.000_kG_RefName\""
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 150,387,"Load Scans:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,24,750,24
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,492,375,492
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 102,514,"Calculate Transmission(s):"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,133,375,133
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 25,147,"Example Scan Name Syntax For Running Code"
	SetDrawEnv fsize= 14,fstyle= 1,textxjust= 2,textyjust= 1
	DrawText 656,38,"Enter Sample Parameters:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 374,24,374,596
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 374,427,750,427
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 474,546,"Make Graph / Add Offset:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,80,375,80
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 485,209,"Calculate Conductivity:"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 446,295,"Calculate Magnetic Susceptibility:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 376,273,749,273
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 376,342,752,342
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 420,369,"Faraday Rotation / Constant Field Folders"
	Button LoadRefScans,pos={7,400},size={180,40},proc=RefsInPhase,title="Load Reference Scans\r In Phase"
	Button LoadRefScans,fSize=12,fStyle=1
	Button LoadSampScans,pos={7,445},size={180,40},proc=SampsInPhase,title="Load Sample Scans\r In Phase"
	Button LoadSampScans,fSize=12,fStyle=1
	Button CalcTrans,pos={7,526},size={180,40},proc=transcalc,title="Calculate Transmission(s)"
	Button CalcTrans,fSize=12,fStyle=1
	Button button4,pos={382,219},size={180,40},proc=ConductivityCalc,title="Calculate Conductivity \r/ Conductance"
	Button button4,fSize=12,fStyle=1
	SetVariable SampleName,pos={7,275},size={204,20},title="\\Z14\\F'Arial'Sample Name?     "
	SetVariable SampleName,value= root:Packages:NickLTHzExtras:Sample_Name
	SetVariable RefName,pos={7,295},size={204,20},title="\\Z14\\F'Arial'Reference Name? "
	SetVariable RefName,value= root:Packages:NickLTHzExtras:Ref_Name
	SetVariable SampleThickness,pos={400,63},size={200,20},title="\\Z14\\F'Arial'Sample Thickness (mm)"
	SetVariable SampleThickness,value= root:Packages:NickLTHzExtras:dsc
	Button button5,pos={565,463},size={180,40},proc=FrequencyCuts,title="Perform Frequency Cuts"
	Button button5,fSize=12,fStyle=1
	SetVariable StepSize,pos={380,498},size={185,20},title="\\Z14\\F'Arial'Step Size (THz)      "
	SetVariable StepSize,value= root:Packages:NickLTHzExtras:Cut_Resolution
	SetVariable StartFreq,pos={380,452},size={185,20},title="\\Z14\\F'Arial'Starting Freq. (THz)"
	SetVariable StartFreq,value= root:Packages:NickLTHzExtras:Freq_Low_Limit
	SetVariable EndFreq,pos={380,475},size={185,20},title="\\F'Arial'\\Z14Ending Freq. (THz) "
	SetVariable EndFreq,value= root:Packages:NickLTHzExtras:Freq_High_Limit
	Button GetFolderInfo,pos={97,324},size={180,40},proc=FolderGrab,title="Get Folder Info"
	Button GetFolderInfo,fSize=12,fStyle=1
	SetVariable SubIndex,pos={400,101},size={200,20},title="\\F'Arial'\\Z14Sub Index of Refraction "
	SetVariable SubIndex,value= root:Packages:NickLTHzExtras:ns
	SetVariable DeltaL,pos={615,101},size={133,20},title="\\Z14\\F'Arial'\\F'Symbol'D\\F'Arial'L (microns)\\Z16"
	SetVariable DeltaL,value= root:Packages:NickLTHzExtras:DeltaL
	PopupMenu SampleType,pos={7,109},size={251,21},proc=MeasurementType,title="\\F'Arial'\\Z14Sample:"
	PopupMenu SampleType,mode=3,popvalue="Thin Film With Substrate Reference",value= #"\"Single Crystal With Aperture Reference;Single Crytal on Substrate With Substrate Reference;Thin Film With Substrate Reference;\""
	SetVariable SubIndex2,pos={400,163},size={200,20},title="\\F'Arial'\\Z14Sub Index of Refraction "
	SetVariable SubIndex2,value= root:Packages:NickLTHzExtras:ns
	SetVariable SampThickness,pos={400,142},size={200,20},title="\\Z14\\F'Arial'Sample Thickness (mm)"
	SetVariable SampThickness,value= root:Packages:NickLTHzExtras:dsc
	SetVariable OffsetAmount,pos={584,562},size={153,20},title="\\Z14\\F'Arial'Offset Amount"
	SetVariable OffsetAmount,value= root:Packages:NickLTHzExtras:OSA
	CheckBox DoGraphs_Checkbox,pos={125,573},size={129,16},proc=CheckGraphs,title="\\Z14\\F'Arial'Generate Graphs"
	CheckBox DoGraphs_Checkbox,value= 0
	CheckBox BLC_Checkbox,pos={8,56},size={100,16},proc=BLC,title="\\Z14\\F'Arial'Temperature"
	CheckBox BLC_Checkbox,value= 0
	CheckBox DoGraphs_Checkbox2,pos={-69,-100},size={129,16},proc=CheckGraphs,title="\\Z14\\F'Arial'Generate Graphs"
	CheckBox DoGraphs_Checkbox2,value= 1
	CheckBox Magnet_Checkbox,pos={139,56},size={110,16},proc=Magnet,title="\\Z14\\F'Arial'Magnetic Field"
	CheckBox Magnet_Checkbox,value= 0
	CheckBox Rotator_Checkbox,pos={273,56},size={93,16},proc=Rotator,title="\\Z14\\F'Arial'Use Rotator"
	CheckBox Rotator_Checkbox,value= 0
	Button LoadSampScans1,pos={187,445},size={180,40},proc=SampsOutPhase,title="Load Sample Scans\r Out Of Phase"
	Button LoadSampScans1,fSize=12,fStyle=1
	Button LoadRefScans1,pos={187,400},size={180,40},proc=RefsOutofPhase,title="Load Reference Scans\r Out of Phase"
	Button LoadRefScans1,fSize=12,fStyle=1
	Button button7,pos={382,297},size={180,40},proc=SuceptibilityCalc,title="Calculate Magnetic \rSusceptibility"
	Button button7,fSize=12,fStyle=1
	SetVariable RefTemp,pos={582,306},size={157,20},title="\\F'Arial'\\Z14Reference Temp"
	SetVariable RefTemp,value= root:Packages:NickLTHzExtras:RefTemp
	Button button8,pos={382,379},size={180,40},proc=FaradayCalc,title="Calculate Faraday Rotation"
	Button button8,fSize=12,fStyle=1
	Button button9,pos={565,379},size={180,40},proc=ConstFieldWaves,title="Generate Constant\rField Folders"
	Button button9,fSize=12,fStyle=1
	Button CalcTrans1,pos={187,526},size={180,40},proc=ConvertToCircBasis,title="Swtich To Circular Basis"
	Button CalcTrans1,fSize=12,fStyle=1
	CheckBox Subn_Checkbox,pos={567,210},size={182,32},proc=Subn,title="\\Z14\\F'Arial'Use Calculated Substrate \rIndex of Refraction"
	CheckBox Subn_Checkbox,value= 0
	SetVariable SubName,pos={584,247},size={156,20},title="\\Z14\\F'Arial'Sub Name? "
	SetVariable SubName,value= root:Packages:NickLTHzExtras:Sub_Name
	Button button0,pos={382,553},size={180,40},proc=MakeAGraph,title="Make Graph"
	Button button0,fSize=12,fStyle=1
EndMacro

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function gets the information of the folder where all the data will be loaded from as well as generates all the data folders in Igor we'll need for our analysis.  Works for both the BLC and Magnet, with and without the rotator.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function GetFolder()

	//---Lets call our global variables---//
		
	//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
	//---This is the case that the "BLC" option is selected on the front panel---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		if (V_Value == 1)
		
			//---Check to see if we are using the rotator in this experiment---//
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			
				//---If we are then create folders for the in phase and out of phase data---//
				if (V_Value == 1)
					
					NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis"
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":'In_Phase'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Time Traces'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'FFTs'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Time Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'FFT Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Transmissions'	
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Mag_Transmissions'	
					
					NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Time Traces'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'FFTs'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Time Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'FFT Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Transmissions'	
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Mag_Transmissions'	
				
				//---If not then we don't need both folders and we'lll just use the in phase ones---//
				else
					
					NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis"
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":'In_Phase'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Time Traces'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'FFTs'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Time Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'FFT Averages'
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Transmissions'	
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Mag_Transmissions'
					
				endif
		endif
		
		//---This is the case that we're using the magnet---//
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		if (V_Value == 1)
			SetDataFolder Root:
			Wave/T Temps
			variable i, j=0, NumTemps;
			
			NumTemps = NumPnts(Temps)
			
			for(i=0; i<NumTemps; i+=1)
					
				//---Checks to see if we're using the rotator---//
				ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					//---This we are, makes in phase and out of phase folders---//
					if (V_Value == 1)
	
						//-----Create the data folders, each temperature gets one of each-----//
						NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis"
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Time Traces'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'FFTs'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Transmissions'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Mag_Transmissions'
						
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase':$Temps[i]
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase':$Temps[i]:'Time Traces'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase':$Temps[i]:'FFTs'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase':$Temps[i]:'Transmissions'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Out_of_Phase':$Temps[i]:'Mag_Transmissions'
						
					//---If not then we don't need both folders---//
					else
					
					//-----Create the data folders, each temperature gets one of each-----//
						NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis"
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Time Traces'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'FFTs'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Transmissions'
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'In_Phase':$Temps[i]:'Mag_Transmissions'
					
					endif
			endfor
		endif

	//---Grabs the file path for the file that was selected---///
	String Path
	GetFileFolderInfo/Q
	Path = S_Path

	//-----Make a string that just has the folder name-----//
	variable folder_ending = strsearch(path,":",Inf,1)
	String folder = path[0,folder_ending]
	
	//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
	base_folder = folder
	
	//---Reports a problem if their is an error or if the user cancels the selection.  Aborts the code---//
	if (V_Flag == -1)
		print "User Cancelled Selection"
		abort
	endif
	
	if (V_Flag != 0 && V_Flag != -1)
		print "File Does Not Exist"
		abort
	endif	
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This is a function I wrote to load all the scans for any experiment.  Folder1 can be either "In_Phase" or "Out_of_Phase".  ScanName is set to either the reference or sample name depending on if you enter in "Ref_Name" 
//  or "Sample_Name" when calling the function, no onther options are valid. Endname is set to either "InP" or "OP" depending on if Folder1 is "In_Phase" or "Out_of_Phase"
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function LoadScans(Folder1, ScanName)

//---Define our parameters---//
	String Folder1, ScanName;
	
//-----String where you enter in the name of the sample and aperature-----//
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	
//---This if statement sets ScanName to be either the reference or the sample name depending on how we're using it---//
	If (cmpstr(ScanName, "Ref_Name") == 0)
		ScanName = Ref_Name
	elseIf (cmpstr(ScanName, "Sample_Name") == 0)
		ScanName = Sample_Name
	else
		abort
	endif

//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----Dummy variables-----//
	String fname, fname2, fname3, fname4;
	
//----Other Strings we'll need----//
	String Folder2, EndName;
	
//---This if statement sets the value of folder2 depending on what folder 1 is, this is so we can load in phase scans in the out of scan fft folder etc---//
	If (cmpstr(Folder1, "In_Phase") == 0)
		Folder2 = "Out_Of_Phase"
		EndName = "InP"
	elseIf (cmpstr(Folder1, "Out_of_Phase") == 0)
		Folder2 = "In_Phase"
		EndName = "OP"
	endif

//---This is the case that the "BLC" option is selected on the front panel---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		if (V_Value == 1)
						
			//---Sets the data folder to the Analysis folder---//
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
			
			//---recognizes the waves that we'll need while running this function---//
			Wave/T Temps
			Wave Num_Ref_scans
			Wave Num_Samp_scans
			Wave NumIter
			Variable NumTemps = numpnts(Temps)
			
			//---This if statement sets the loop iteration to be Num_Ref_Scans or Num_Samp_Scans based on "ScanName"---//
			If (cmpstr(ScanName, Ref_Name) == 0)
				Duplicate/O Num_Ref_Scans, NumIter
			elseIf (cmpstr(ScanName, Sample_Name) == 0)
				Duplicate/O Num_Samp_Scans, NumIter
			else
				abort
			endif
		
			//---Loop that will load all the scans---//
			Variable i, j;
			for(i = 0; i < NumTemps; i += 1)
				for( j = 1; j <= NumIter[i]; j += 1)
				
					//---Load all the scans---//
					fname = (ScanName + "_" + Temps[i] + "_" + num2str(j))
					SingleWaveLoad(base_folder,fname)
					fft_gen(fname)
		
					//----Create the general wave name for this temperature----//
					fname = ScanName + "_" + Temps[i]
					fname2 = ScanName + "TAvg_" + EndName
				
				//-----Add all sample scans together-----//
					//---Set the jth version of this temperature as the current wave to add to the average---//
					Wave Current_ref = $(fname + "_" + num2str(j))
					
					//----Create the averaged scans if first sample scans, otherwise add the scan to the average----//
					if(j==1)
						Duplicate/O $(fname + "_" + num2str(j)), $(fname + "_Avg_" + EndName)
						Wave ref_TAvg = $(fname + "_Avg_" + EndName)
					else				
						ref_TAvg += Current_ref
					endif
				endfor
					//-----Divide by the number of scans to get the averaged pulse-----//
					ref_TAvg /= NumIter[i]
					
					//-----Do the FFT of the pulse and then do some wave clean up---//
					fft_gen(fname + "_Avg_" + EndName)
					
				endfor
				
			for(i = 0; i < NumTemps; i += 1)
				for( j = 1; j <= NumIter[i]; j += 1)
						
						fname = (ScanName + "_" + Temps[i] + "_" + num2str(j))
						Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'Time Traces':$fname
						KillWaves $fname
						
						fname = (ScanName + "_" + Temps[i] + "_" + num2str(j) + "_FFT")
						Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):FFTs:$(fname)
						KillWaves $fname
				
				endfor
						fname = (ScanName + "_" + Temps[i] + "_Avg_" + EndName + "_FFT")
						Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'FFT Averages':$fname
						
						ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
						if (V_Value==1)
							Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder2):'FFT Averages':$fname
						endif
						KillWaves $fname
						
						fname = (ScanName + "_" + Temps[i] + "_Avg_" + EndName)
						Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'Time Averages':$fname
						KillWaves $fname
			endfor
			
			KillWaves/Z NumIter
		endif		
			
//---This is the case that we're using the magnet---//
	ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
	if (V_Value == 1)
		//---Define some variables---//
		variable NumFields;
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis":
			
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
	
				fname = (Temps[i] + "_" + Fields[j] + "0_kG_" + ScanName + ".txt")
				fname2 = (Temps[i] + "_" + Fields[j] + "0_kG_" + ScanName)
				fname3 = (ScanName + "_" + Temps[i] + "_" + Fields[j] + "kG_" + EndName)
				fname4 = (fname3 + "_FFT")			
									
				// Loads the scans based on the parameters in Chris's panel.
				SingleWaveLoad(base_folder, fname)
				
				// Renames the scans
				Duplicate/O $fname2, $fname3
				Killwaves/Z $fname2
				
				// Now take the FFT of the scan
				fft_gen(fname3)
					
				//---Duplicate the Time Trace and move it to the time traces folder---//
				Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:'Time Traces':$fname3
				KillWaves/Z $fname3
				
				// Duplicate/O the FFT and move it to the in Phase FFT folder and the Out of Phase FFT Folder
				Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:'FFTs':$fname4
				
				//---Checks to see if we're using the rotator---//
				ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)
						Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":$(Folder2):$Temps[i]:'FFTs':$fname4
					endif
					
				// Delete the FFT From the time traces folder.	
				KillWaves/Z $fname4
			endfor	
		endfor
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This is a function I wrote to calculate the transmission for either the BLC or the magnet, with or without the rotator.  Folder1 can be either "In_Phase" or "Out_Of_Phase".
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function CalcTrans(Folder1)
//----Define our function parameters---//
	String Folder1;

//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
	
//-----Dummy variables-----//
	String fname, fname2, fname3, fname4, fname5, fname6, fname7, fname8;
	
//---Some other strings we'll need for the calculation---//
	String Direction, EndName1, EndName2;
	
//---This if statement sets the value of folder2 depending on what folder 1 is, this is so we can load in phase scans in the out of scan fft folder etc---//
	If (cmpstr(Folder1, "In_Phase") == 0)
		Direction = "xx"
		EndName1 = "InP"
		EndName2 = "OP"
	elseIf (cmpstr(Folder1, "Out_of_Phase") == 0)
		Direction = "yx"
		EndName1 = "OP"
		EndName2 = "InP"
	endif
	
//---This is the case that the "BLC" option is selected on the front panel---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		if (V_Value == 1)
						
			//---Sets the data folder to the Analysis folder---//
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
			
			//---recognizes the waves that we'll need while running this function---//
			Wave/T Temps
			Variable NumTemps = numpnts(Temps)
			Variable delta, i, j
			
			for(i = 0; i < NumTemps; i+=1)
				
				//------Calculate the transmission-----//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'FFT Averages'
				Duplicate/O $(Sample_Name + "_" + Temps[i] + "_Avg_" + EndName1 + "_FFT"), $("T" + Direction + "_" + Sample_Name + "_" + Temps[i])
				Wave/C top = $("T" + Direction + "_" + Sample_Name + "_" + Temps[i])
				Wave/C bottom = $(Ref_Name + "_" + Temps[i] + "_Avg_InP_FFT")
				top /= bottom
				
				//---Duplicate them into the correct folders and then delete them---// 
				fname = ("T" + Direction + "_" + Sample_Name + "_" + Temps[i])
				fname2 = ("Mag_T" + Direction + "_" + Sample_Name + "_" + Temps[i])
				Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'Transmissions':$fname
				Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'Mag_Transmissions':$fname2
				KillWaves/Z $fname
				
				//---Go to the Mag_Transmissions folder and convert the waves to polar form then redimension them as real---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):'Mag_Transmissions'
				Wave/C Fake_wave = $fname2
				Fake_wave = r2polar(Fake_wave)
				Redimension/R Fake_wave
					
			endfor
		endif
		
//---This is the case that we're using the magnet---//
	ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		if (V_Value == 1)
		
		//---Define some variables---//
		variable NumFields;
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis":
			
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)

				//---Lets first compute txx, so set the folder to be the in phase FFTs---//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:'FFTs'
			
				//---First strings are the in phase and out of phase sample scans---//
				fname = (Sample_Name + "_" + Temps[i] + "_" + Fields[j] + "kG_InP_FFT")	
				fname2 = (Sample_Name + "_" + Temps[i] + "_" + Fields[j] + "kG_OP_FFT")	
				
				//---Next two are the in phase and our of phase aperture scans---//
				fname3 =  (Ref_Name + "_" + Temps[i] + "_" + Fields[j] + "kG_InP_FFT")	
				fname4 =  (Ref_Name + "_" + Temps[i] + "_" + Fields[j] + "kG_OP_FFT")	
				
				//---Then the transmission and magnitude of the transmission---//
				fname5 = ("T" + Direction + "_" + Temps[i] + "_" + Fields[j] + "kG")
				fname6 = ("Mag_"+fname5)
				
				//---Duplicate one of the FFTs and call it the tranmission so we have that wave---//
				Duplicate/O $fname, $fname5
				
				//---Name some dummy waves to hold the waves for our txx and tyx calculation---//
				Wave/C E_OutX 	= $fname
				Wave/C E_OutY 	= $fname2
				Wave/C E_InX 	= $fname3
				Wave/C E_InY 	= $fname4
				Wave/C T		= $fname5
					
				if (cmpstr(Direction, "xx") == 0)
					T = E_OutX / E_InX
					//T = (E_OutX*E_InX + E_OutY*E_InY) / (E_InX^2 + E_InY^2)
				elseif (cmpstr(Direction, "yx") == 0)
					//T = E_OutY / E_InX
					T = (E_OutY*E_InX - E_OutX*E_InY) / (E_InX^2 + E_InY^2)							
				endif
				
				//---Now Duplicate/O the transmission and store it in the transmissions folder---//
				Duplicate/O $fname5, root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:Transmissions:$fname5
				
				//---Do the same thing for the magnitude of the transmissions---//
				Duplicate/O $fname5,  root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:Mag_Transmissions:$fname6
				
				//---Detele the original from the FFTs folder---//
				KillWaves/Z $fname5
				
				//---Now go to the Mag Transmissions Folder---//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:Mag_Transmissions
				
				//----New Dummy wave that we'll convert to polar and redimension so that the magnitude is a real wave---//
				Wave/C Fake_wave = $fname6
				Fake_wave = r2polar(Fake_wave)
				Redimension/R Fake_wave
			endfor
		endfor
	endif
End	
	
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function takes the transmission data that we have already computed and turns it into the circular basis by a simple linear transformation.  Only works when using the rotator since we need Tyx for the transformation.
// Really will only be useful in zero field with the rotator or on the magnet in Faraday geometry with the rotator.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
Function CircBasisTrans()

//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----Dummy variables-----//
	String fname, fname2, fname3, fname4, fname5, fname6
	
	variable i,k,j;
	
	//---Function only runs when the rotator is being used---//
	ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
	if (V_Value ==1)
	
		//---This is the case that the "BLC" option is selected on the front panel---//
			ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
			if (V_Value == 1)
						
				//---Sets the data folder to the Analysis folder---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
				
				//---recognizes the waves that we'll need while running this function---//
				Wave/T Temps
				Variable NumTemps = numpnts(Temps)
				
				//---Creates folders for the left and right circular basis---//
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Mag_Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Mag_Transmissions
		
				
				//---Loop that runs over our temperatures---//
				for(i=0; i<NumTemps; i+=1)
				
					//---These strings hold the linear basis transmissions---//
					fname = ("Txx_" + Sample_Name + "_" + Temps[i])
					fname2 = ("Tyx_" + Sample_Name + "_" + Temps[i])
					
					//----These strings hold the circular basis transmissions---//
					fname3 = ("Tr_" + Sample_Name + "_" + Temps[i])
					fname4 = ("Tl_" + Sample_Name + "_" + Temps[i])
					
					//----These strings hold the magnitude of the circular basis transmissions---//
					fname5 = ("Mag_Tr_" + Sample_Name + "_" + Temps[i])
					fname6 = ("Mag_Tl_" + Sample_Name + "_" + Temps[i])

					//----First go to the in phase transmissions and Duplicate Txx---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Transmissions
					Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:$fname
					
					//---Then go to the our of phase transmissions and Duplicate Tyx---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Transmissions
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:$fname2
					
					//---Duplicate some wave and rename them as the ciruclar basis transmissions---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:$fname3
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:$fname4
					
					//---Declare some dummy waves for the calculation---//
					Wave/C DummyWave = $fname
					Wave/C DummyWave2 = $fname2
					Wave/C DummyWave3 = $fname3
					Wave/C DummyWave4 = $fname4
					
					//---Convert to the circular basis---//
					DummyWave3 = DummyWave + sqrt(-1)*Dummywave2
					DummyWave4 = DummyWave - sqrt(-1)*Dummywave2
					
					//---Duplicate/O the circular transmissions and store them in the correct folder---//
					Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Transmissions:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$fname4
					
					//---Duplicate/O the circular transmissions and store them in the magnitude Folder as well---//
					Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Mag_Transmissions:$fname5
					Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Mag_Transmissions:$fname6

					//---Delete the linear transmissions---//
					KillWaves/Z $fname, $fname2, $fname3, $fname4
					
					//---Now lets go to the magnitude folders and convert those complex transmissions to real waves---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Mag_Transmissions
					Wave/C Fake_wave = $fname5
					Fake_wave = r2polar(Fake_wave)
					Redimension/R Fake_wave
					
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Mag_Transmissions
					Wave/C Fake_wave = $fname6
					Fake_wave = r2polar(Fake_wave)
					Redimension/R Fake_wave
				endfor							
			endif
			
		//---This is the case that we're using the magnet---//
			ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
			if (V_Value == 1)
			
			//---Creates a general folder for the circular basis data---//
			NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis
			
			//---Define some variables---//
			variable NumFields;
	
			SetDataFolder root:$Sample_Name + "_InField_Analysis":
			Wave/T Temps
			NumTemps = NumPnts(Temps)
		
			for (i=0; i<NumTemps; i+=1)
		
				//---Set the data folder to the main folder to find "Fields_Temps" waves
				SetDataFolder root:$Sample_Name + "_InField_Analysis":
				
				fname = "Fields_" + Temps[i]
				Wave/T Fields = $fname
				NumFields = NumPnts($fname)
				
				//---Creates the folders we'll need for converting to the circular basis---//
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Mag_Transmissions
				NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Mag_Transmissions
				
			
				//--- j runs over the number of Fields we have---//
				for (j=0; j<NumFields; j+=1)
				
					//---These strings hold the linear basis transmissions---//
					fname = ("Txx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname2 = ("Tyx_" + Temps[i] + "_" + Fields[j] + "kG")
					
					//----These strings hold the circular basis transmissions---//
					fname3 = ("Tr_" + Temps[i] + "_" + Fields[j] + "kG")
					fname4 = ("Tl_" + Temps[i] + "_" + Fields[j] + "kG")
					
					//---These strings hold the magnitude of the circular basis transmissions---//
					fname5 = ("Mag_Tr_" + Temps[i] + "_" + Fields[j] + "kG")
					fname6 = ("Mag_Tl_" + Temps[i] + "_" + Fields[j] + "kG")

					//---First go to the in phase transmissions and Duplicate/O Txx---//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions
					Duplicate/O $fname, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:$fname
					
					//---Then go to the our of phase transmissions and Duplicate/O Tyx---//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Transmissions
					Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:$fname2
					
					//---Duplicate/O some wave and rename them as the ciruclar basis transmissions---//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]
					Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:$fname3
					Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:$fname4
					
					//---Declare some dummy waves for the calculation---//
					Wave/C DummyWave = $fname
					Wave/C DummyWave2 = $fname2
					Wave/C DummyWave3 = $fname3
					Wave/C DummyWave4 = $fname4
					
					//---Convert to the circular basis---//
					DummyWave3 = DummyWave + sqrt(-1)*Dummywave2
					DummyWave4 = DummyWave - sqrt(-1)*Dummywave2
					
					//---Duplicate/O the circular transmissions and store them in the correct folder---//
					Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$fname4
					
					//---Duplicate/O the circular transmissions and store them in the magnitude Folder as well---//
					Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Mag_Transmissions:$fname5
					Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Mag_Transmissions:$fname6
					
					//---Delete the linear transmissions---//
					KillWaves/Z $fname, $fname2, $fname3, $fname4
					
					//---Now lets go to the magnitude folders and convert those complex transmissions to real waves---//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Mag_Transmissions
					Wave/C Fake_wave = $fname5
					Fake_wave = r2polar(Fake_wave)
					Redimension/R Fake_wave
					
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Mag_Transmissions
					Wave/C Fake_wave = $fname6
					Fake_wave = r2polar(Fake_wave)
					Redimension/R Fake_wave
				
				endfor	
			endfor		
		endif
		
		else
			Print "You must use the rotator to convert to the circular basis!"
		endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the conductivity for a single crystal for either the BLC or Magnet, with or without the rotator.  Conductivity for a thin film sample will be written separately.  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function SCCond()
//-----Folder where data is loaded from-----//
SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----String where you enter in the name of the sample and aperature-----//
SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name

//-----Thickness of the sample in mm's------//
NVAR dsc = root:Packages:NickLTHzExtras:dsc

//-----PopVal sets the type of sample being measured so we know which version of NR code to run------//
NVAR PopVal = root:Packages:NickLTHzExtras:PopVal

//-----Substrate Index of Refraction, we'll need this for a single crystal on a substrate with a substrate reference------//
NVAR ns = root:Packages:NickLTHzExtras:ns

//---Checks to see if we're running the BLC---//
ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
	if (V_Value == 1)
	
		//---Sets the data folder to the Analysis folder just in case---//
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
		
		//---Recognizes the wave that we'll need while running this function---//
		Wave/T Temps
		Wave num_samp_scans
		Variable NumTemps = numpnts(Temps)
		
		//-----Dummy variables-----//
		String fname, fname2;
		
		//---Creates folders for the phase and the real and imaginary part of sigma and n for the in phase and out of phase folders---//
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:k
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Sigma1
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Sigma2
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Phase
		
		//------Calculate the index of refraction of the sample------//
		variable delta, i, j
		
		//--The number of Netwon-Rhapson iterations---//
		Variable num_iter = 10 
		
		//---Loop runs over the temperatures in the experiment---//
		for(i = 0; i < NumTemps; i+=1)
		
			//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Transmissions'
			index_calc("Txx_" + Sample_Name + "_" + Temps[i])
			Wave transmission, phase, freq
			
			//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
			CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
			Wave W_coef
			phase -= W_coef[0]
			
			//-----Run the newton-raphson code to find the index-----//
			delta = DimDelta($("Txx_" + Sample_Name + "_" + Temps[i]),0)
			
			//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"---//
			//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
			
			if (PopVal==1)
				//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
				TDNR(transmission, phase, freq, dsc, num_iter)
			elseif (PopVal==2)
				//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
				TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
			endif
			
			//---Declare n and k and then set their scale---//
			Wave n,k
			SetScale/P x 0,(delta),"", n,k
			
			//-----Calculate the conductivity-----//
			Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
			Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
			sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
			
			//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
			Duplicate/O $"n",  root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n:$("nxx_" + Sample_Name + "_" + Temps[i])
			Duplicate/O $"k", root:$Sample_Name + "_ZeroField_Analysis":In_Phase:k:$("kxx_" + Sample_Name + "_" + Temps[i])
			Duplicate/O $"phase", root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Phase:$("Thetaxx_" + Sample_Name + "_" + Temps[i])
			Duplicate/O $"sigma1", root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Sigma1:$("sigma1xx_" + Sample_Name + "_" + Temps[i])
			Duplicate/O $"sigma2", root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Sigma2:$("sigma2xx_" + Sample_Name + "_" + Temps[i])
			
		endfor
		
		//-----A bit of wave clean up-----//
		Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
		KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
		Wave sigma1, sigma2
		KillWaves sigma1, sigma2
			
		ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
		if (V_Value == 1)
		
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:k
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Sigma1
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Sigma2
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Phase
			
			//---Now repeat the entire process for the out of phase waves---//
			for(i=0; i<NumTemps; i+=1)
			
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Transmissions'
				index_calc("Tyx_" + Sample_Name + "_" + Temps[i])
				Wave transmission, phase, freq
				
				//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
				CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
				Wave W_coef
				phase -= W_coef[0]
				
				//-----Run the Newton-Raphson code to find the index-----//
				delta = DimDelta($("Tyx_" + Sample_Name + "_" + Temps[i]),0)
				
				//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
				//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
				
				if (PopVal==1)
					//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
					TDNR(transmission, phase, freq, dsc, num_iter)
				elseif (PopVal==2)
					//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
					TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
				endif
				
				//---Declare n and k and then set their scale---//
				Wave n,k
				SetScale/P x 0,(delta),"", n,k
				
				//-----Calculate the conductivity-----//
				Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
				Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
				sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
				
				//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
				Duplicate/O $"n",  root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n:$("nyx_" + Sample_Name + "_" + Temps[i])
				Duplicate/O $"k", root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:k:$("kyx_" + Sample_Name + "_" + Temps[i])
				Duplicate/O $"phase", root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Phase:$("Thetayx_" + Sample_Name + "_" + Temps[i])
				Duplicate/O $"sigma1", root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Sigma1:$("sigma1yx_" + Sample_Name + "_" + Temps[i])
				Duplicate/O $"sigma2", root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Sigma2:$("sigma2yx_" + Sample_Name + "_" + Temps[i])
				
			endfor
		
			//-----A bit of wave clean up-----//
			Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			Wave sigma1, sigma2
			KillWaves sigma1, sigma2
		endif
	endif
	
	//---Now the case that we're running the magnet with or without the rotator---//
	ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
	if (V_Value == 1)
		
		//---Define some variables---//
		variable NumFields;
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis":
			
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
			
			//---Creates folders for the complex conductivity, complex n, and phase---//
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:k
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Sigma1
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Sigma2
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Phase
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
		
				SetDataFolder root:$Sample_Name + "_InField_Analysis":
				
				//------Calculate the index of refraction of the sample------//
				//--The number of Netwon-Rhapson iterations---//
				 num_iter = 10 
			
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions
				index_calc("Txx_" + Temps[i] + "_" + Fields[j] + "kG")
				Wave transmission, phase, freq
				
				//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
				CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
				Wave W_coef
				phase -= W_coef[0]
				
				//-----Run the newton-raphson code to find the index-----//
				delta = DimDelta($("Txx_" + Temps[i] + "_" + Fields[j] + "kG"),0)
				
				//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
				//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
				
				if (PopVal==1)
					//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
					TDNR(transmission, phase, freq, dsc, num_iter)				
				elseif (PopVal==2)
					//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
					TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)							
				endif
				
				//---Declare n and k and then set their scale---//
				Wave n,k
				SetScale/P x 0,(delta),"", n,k
				
				//-----Calculate the conductivity-----//
				Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
				Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
				sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
				
				//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
				Duplicate/O $"n", root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n:$("nxx_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $"k",  root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:k:$("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $"phase", root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Phase:$("Thetaxx_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $"sigma1",  root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Sigma1:$("sigma1xx_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $"sigma2", root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Sigma2:$("sigma2xx_" + Temps[i] + "_" + Fields[j] + "kG")
					
			endfor
			
			//-----A bit of wave clean up-----//
			Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			Wave sigma1, sigma2
			KillWaves sigma1, sigma2
			
		ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
		if (V_Value == 1)
		
			//---Creates folders for the complex conductivity, complex n, and phase---//
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:n
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:k
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Sigma1
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Sigma2
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Phase
				
			
			//---Do the same thing for the out of phase waves---//
			for(j = 0; j < NumFields; j+=1)
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Transmissions
					index_calc("Tyx_" + Temps[i] + "_" + Fields[j] + "kG")
					Wave transmission, phase, freq
					
					//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
					CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
					Wave W_coef
					phase -= W_coef[0]
					
					//-----Run the newton-raphson code to find the index-----//
					delta = DimDelta($("Tyx_" + Temps[i] + "_" + Fields[j] + "kG"),0)
					
					//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
					//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
					
					if (PopVal==1)
						//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
						TDNR(transmission, phase, freq, dsc, num_iter)
					elseif (PopVal==2)
						//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
						TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
					endif
					
					//---Declare n and k and then set their scale---//
					Wave n,k
					SetScale/P x 0,(delta),"", n,k
					
					//-----Calculate the conductivity-----//
					Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
					Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
					sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
					
					//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
					Duplicate/O $"n", root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n:$("nyx_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $"k",  root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:k:$("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $"phase", root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Phase:$("Thetayx_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $"sigma1",  root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Sigma1:$("sigma1yx_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $"sigma2", root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Sigma2:$("sigma2yx_" + Temps[i] + "_" + Fields[j] + "kG")
					
			endfor
			
				//-----A bit of wave clean up-----//
				Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
				KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
				Wave sigma1, sigma2
				KillWaves sigma1, sigma2
		endif
	endfor
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the conductivity for a single crystal for either the BLC or Magnet, in the circular basis.  Only works if we're using the roator and if the Circular Basis data folder exists.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function CircBasisSCCond()
//-----Folder where data is loaded from-----//
SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----String where you enter in the name of the sample and aperature-----//
SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name

//-----Thickness of the sample in mm's------//
NVAR dsc = root:Packages:NickLTHzExtras:dsc

//-----PopVal sets the type of sample being measured so we know which version of NR code to run------//
NVAR PopVal = root:Packages:NickLTHzExtras:PopVal

//-----Substrate Index of Refraction, we'll need this for a single crystal on a substrate with a substrate reference------//
NVAR ns = root:Packages:NickLTHzExtras:ns
	
	//---Checks to see if we've already created the circular basis folder.
	if (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:Circular_Basis") == 1)
		
		//---Checks to see if we're running the BLC---//
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
			if (V_Value == 1)
			
				//---Sets the data folder to the Analysis folder just in case---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
				
				//---Recognizes the wave that we'll need while running this function---//
				Wave/T Temps
				Wave num_samp_scans
				Variable NumTemps = numpnts(Temps)
				
				//-----Dummy variables-----//
				String fname, fname2;
				
				//---Creates folders for the Right Hand circular Data---//
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:k
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Sigma1
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Sigma2
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Phase

				
				//------Calculate the index of refraction of the sample------//
				variable delta, i, j
				
				//--The number of Netwon-Rhapson iterations---//
				Variable num_iter = 10 
				
				//---Loop runs over the temperatures in the experiment---//
				for(i = 0; i < NumTemps; i+=1)
				
					//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:'Transmissions'
					index_calc("Tr_" + Sample_Name + "_" + Temps[i])
					Wave transmission, phase, freq
					
					//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
					CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
					Wave W_coef
					phase -= W_coef[0]
					
					//-----Run the newton-raphson code to find the index-----//
					delta = DimDelta($("Tr_" + Sample_Name + "_" + Temps[i]),0)
					
					//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"---//
					//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
					
					if (PopVal==1)
						//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
						TDNR(transmission, phase, freq, dsc, num_iter)
					elseif (PopVal==2)
						//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
						TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
					endif
					
					//---Declare n and k and then set their scale---//
					Wave n,k
					SetScale/P x 0,(delta),"", n,k
					
					//-----Calculate the conductivity-----//
					Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
					Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
					sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
					
					//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
					Duplicate/O $"n",  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n:$("nr_" + Sample_Name + "_" + Temps[i])
					Duplicate/O $"k", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:k:$("kr_" + Sample_Name + "_" + Temps[i])
					Duplicate/O $"phase", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Phase:$("Thetar_" + Sample_Name + "_" + Temps[i])
					Duplicate/O $"sigma1", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Sigma1:$("sigma1r_" + Sample_Name + "_" + Temps[i])
					Duplicate/O $"sigma2", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Sigma2:$("sigma2r_" + Sample_Name + "_" + Temps[i])
					
				endfor
				
				//-----A bit of wave clean up-----//
				Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
				KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
				Wave sigma1, sigma2
				KillWaves sigma1, sigma2

				//---Now do everything again for the left hand data---//
				//---Creates folders for the  Left Hand Circular Data---//
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:k
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Sigma1
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Sigma2
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Phase
					
					//---Now repeat the entire process for the out of phase waves---//
					for(i=0; i<NumTemps; i+=1)
					
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:'Transmissions'
						index_calc("Tl_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the Newton-Raphson code to find the index-----//
						delta = DimDelta($("Tl_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Calculate the conductivity-----//
						Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
						Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
						sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
						
						//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $"n",  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n:$("nl_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $"k", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:k:$("kl_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $"phase", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Phase:$("Thetal_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $"sigma1", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Sigma1:$("sigma1l_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $"sigma2", root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Sigma2:$("sigma2l_" + Sample_Name + "_" + Temps[i])
						
					endfor
				
					//-----A bit of wave clean up-----//
					Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					Wave sigma1, sigma2
					KillWaves sigma1, sigma2
			endif
	endif
	
	//---Checks to see if we've already created the circular basis folder.
	if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:Circular_Basis") == 1)
			
			//---Now the case that we're running the magnet with or without the rotator---//
			ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
			if (V_Value == 1)
				
				//---Define some variables---//
				variable NumFields;
			
				SetDataFolder root:$Sample_Name + "_InField_Analysis":
				Wave/T Temps
				NumTemps = NumPnts(Temps)
				
				for (i=0; i<NumTemps; i+=1)
				
					//---Set the data folder to the main folder to find "Fields_Temps" waves
					SetDataFolder root:$Sample_Name + "_InField_Analysis":
					
					fname = "Fields_" + Temps[i]
					Wave/T Fields = $fname
					NumFields = NumPnts($fname)
					
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:k
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Sigma1
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Sigma2
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Phase
				
					//--- j runs over the number of Fields we have---//
					for (j=0; j<NumFields; j+=1)
				
						SetDataFolder root:$Sample_Name + "_InField_Analysis":
						
						//------Calculate the index of refraction of the sample------//
						//--The number of Netwon-Rhapson iterations---//
						 num_iter = 10 
					
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions
						index_calc("Tr_" + Temps[i] + "_" + Fields[j] + "kG")
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("Tr_" + Temps[i] + "_" + Fields[j] + "kG"),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)				
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)							
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Calculate the conductivity-----//
						Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
						Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
						sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
						
						//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $"n", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n:$("nr_" + Temps[i] + "_" + Fields[j] + "kG")
						Duplicate/O $"k",  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:k:$("kr_" + Temps[i] + "_" + Fields[j] + "kG")
						Duplicate/O $"phase", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Phase:$("Thetar_" + Temps[i] + "_" + Fields[j] + "kG")
						Duplicate/O $"sigma1",  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Sigma1:$("sigma1r_" + Temps[i] + "_" + Fields[j] + "kG")
						Duplicate/O $"sigma2", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Sigma2:$("sigma2r_" + Temps[i] + "_" + Fields[j] + "kG")
							
					endfor
					
					//-----A bit of wave clean up-----//
					Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					Wave sigma1, sigma2
					KillWaves sigma1, sigma2
					
				ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
				if (V_Value == 1)
				
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:k
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Sigma1
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Sigma2
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Phase
						
					
					//---Do the same thing for the out of phase waves---//
					for(j = 0; j < NumFields; j+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions
							index_calc("Tl_" + Temps[i] + "_" + Fields[j] + "kG")
							Wave transmission, phase, freq
							
							//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
							CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
							Wave W_coef
							phase -= W_coef[0]
							
							//-----Run the newton-raphson code to find the index-----//
							delta = DimDelta($("Tl_" + Temps[i] + "_" + Fields[j] + "kG"),0)
							
							//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
							//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
							
							if (PopVal==1)
								//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
								TDNR(transmission, phase, freq, dsc, num_iter)
							elseif (PopVal==2)
								//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
								TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
							endif
							
							//---Declare n and k and then set their scale---//
							Wave n,k
							SetScale/P x 0,(delta),"", n,k
							
							//-----Calculate the conductivity-----//
							Duplicate/O $"n", $"sigma1", $"sigma2"; Wave sig1 = $"sigma1"; Wave sig2 = $"sigma2";
							Wave freq = $"freq"; Wave n = $"n"; Wave k = $"k";
							sig1 = freq / 2 * ( k*n ); sig2 = - freq / 2 * ( (n*n - k*k)/2);
							
							//-----Sort n, k, sigma into the appropriate folders labeled by the current temperature-----//
							Duplicate/O $"n", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n:$("nl_" + Temps[i] + "_" + Fields[j] + "kG")
							Duplicate/O $"k",  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:k:$("kl_" + Temps[i] + "_" + Fields[j] + "kG")
							Duplicate/O $"phase", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Phase:$("Thetal_" + Temps[i] + "_" + Fields[j] + "kG")
							Duplicate/O $"sigma1",  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Sigma1:$("sigma1l_" + Temps[i] + "_" + Fields[j] + "kG")
							Duplicate/O $"sigma2", root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Sigma2:$("sigma2l_" + Temps[i] + "_" + Fields[j] + "kG")
							
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						Wave sigma1, sigma2
						KillWaves sigma1, sigma2
				endif
			endfor
			endif
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the conductivity for a thin film for either the BLC or Magnet, with or without the rotator.  The substrate index of refraction is treated as a real constant here, could change it...
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function TFCond()

//---Folder where data is loaded from---//
SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//---String where you enter in the name of the sample, aperature and sub if chosen---//
SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
SVAR Sub_Name = root:Packages:NickLTHzExtras:Sub_Name

//---Substrate Index of Refraction, only used if the "use calculated substrate index of refraction" option is NOT chosen---//
NVAR ns = root:Packages:NickLTHzExtras:ns

//---DeltaL is difference in thicnkess between reference and sample substrate in microns----//
NVAR DeltaL = root:Packages:NickLTHzExtras:DeltaL

//---Dummy variables---//
String fname, fname2, fname3, fname4, fname5, fname6;
Variable i, j;
Variable Zo = 376.7, c = 3

	//---This is the case that we're running the BLC---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
	if (V_Value == 1)
	
		//---Sets the data folder to the Analysis folder just in case---//
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
		
		//---Recognizes the wave that we'll need while running this function---//
		Wave/T Temps
		Wave num_samp_scans
		Variable NumTemps = numpnts(Temps)
		
			//---Creates folders for the complex conductance of the thin film---//
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Conductance1
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Conductance2
				
			for(i = 0; i < NumTemps; i+=1)
				
				//--Strings for the in phase transmissions--//
				fname = "Txx_" + Sample_Name + "_" + Temps[i]
				
				//--Strings for the in phase Conductance--//
				fname2 = "Gxx_" + Sample_Name + "_"+ Temps[i]
				fname3 =  "G1xx_" + Sample_Name + "_"+ Temps[i]
				fname4 =  "G2xx_" + Sample_Name + "_"+ Temps[i]
				
				//---Now go to the transmissions folder to copy the waves to make conductance waves---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Transmissions'
				Duplicate/O $fname, $fname2
				Duplicate/O $fname, $fname3
				Duplicate/O $fname, $fname4
				
				Redimension/R $fname3, $fname4
				
				//---Makes a wave for the frequency--//
				Variable delta = DimDelta($fname,0)
				Make/O/N=(numpnts($fname)), freq
				for(j = 0; j<numpnts(freq); j+=1)
					freq[j] = delta*j
				endfor
				SetScale/P x 0,(delta),"", freq
				
				//---dummy waves for calculations---//
				Wave/C Trans = $fname
				Wave/C Cond = $fname2
				Wave Cond1 = $fname3
				Wave Cond2 = $fname4
				
				//---This is the option that we're using an already calculated index of refraction of the substrate via this analysis palette---//
				//---The analysis folder of the substrate has to exist prior to trying to analyze the thin films.  So first run the entire code on the sub as a single crystal---//
				//---Then run the code on the thin film---//
				ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
				If (V_Value == 1)
					
					//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":In_Phase:n
					Duplicate/O $("nxx_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Transmissions:$("nxx_" + Sub_Name + "_" + Temps[i])
					
					//---Do the same thing with k, put it in the Sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":In_Phase:k
					Duplicate/O $("kxx_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Transmissions:$("kxx_" + Sub_Name + "_" + Temps[i])
					
					//---Now set the folder to the sample's transmissions's folder----//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:'Transmissions'
					
					//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
					Duplicate/O $fname, NSub
					
					//---Combine the n and k of the substrate into a single complex wave called NSub---//
					Wave/C Dummywave = Nsub
					Wave Dummywave1 = $("nxx_" + Sub_Name + "_" + Temps[i])
					Wave Dummywave2 = $("kxx_" + Sub_Name + "_" + Temps[i])
					
					Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
				
					//---Calculate the conplex conductance using the substrates complex index of refraction---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
					
					//---Now some wave clean up---//
					//KillWaves/Z Nsub, DummyWave1, DummyWave2
				else
					//---Calculate the conplex conductance---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
				endif
				
				//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
				Cond1 = Real(Cond)
				Cond2 = Imag(Cond)

				//---Stores them in the appropriate folder---//
				Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Conductance1:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Conductance2:$fname4
				
				//---Some more wave clean up---//
				KillWaves/Z $fname2, $fname3, $fname4, freq
			
			endfor
			
			//---Now repeat the entire process for the out of phase components---//
			//---This is the case that we're using the rotator
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			if (V_Value == 1)
				
				//---Creates folders for the complex conductance of the thin film---//
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Conductance1
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Conductance2
					
				for(i = 0; i < numpnts(Temps); i+=1)
					
					//--Strings for the in phase transmissions--//
					fname = "Tyx_" + Sample_Name + "_" + Temps[i]
					
					//--Strings for the in phase Conductance--//
					fname2 = "Gyx_" + Sample_Name + "_"+ Temps[i]
					fname3 =  "G1yx_" + Sample_Name + "_"+ Temps[i]
					fname4 =  "G2yx_" + Sample_Name + "_"+ Temps[i]
					
					//---Now go to the transmissions folder to copy the waves to make conductance waves---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:'Transmissions'
					Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Transmissions:$fname2
					
					//---Makes a wave for the frequency--//
					delta = DimDelta($fname,0)
					Make/O/N=(numpnts($fname)), freq
					for(j = 0; j<numpnts(freq); j+=1)
						freq[j] = delta*j
					endfor
					SetScale/P x 0,(delta),"", freq
					
					//----dummy waves for calculations---//
					Wave/C Cond = $fname2
					Wave/C Trans = $fname
					
					//---Calculate the Complex Conductance---//
					ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
					If (V_Value == 1)
					
						//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Out_Of_Phase:n
						Duplicate/O $("nyx_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Transmissions:$("nyx_" + Sub_Name + "_" + Temps[i])
						
						//---Do the same thing with k, put it in the Sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Out_Of_Phase:k
						Duplicate/O $("kyx_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Transmissions:$("kyx_" + Sub_Name + "_" + Temps[i])
						
						//---Now set the folder to the sample's transmissions's folder----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:'Transmissions'
						
						//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
						Duplicate/O $fname, NSub
						
						//---Combine the n and k of the substrate into a single complex wave called NSub---//
						Wave/C Dummywave = Nsub
						Wave Dummywave1 = $("nyx_" + Sub_Name + "_" + Temps[i])
						Wave Dummywave2 = $("kyx_" + Sub_Name + "_" + Temps[i])
						
						Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
					
						//---Calculate the conplex conductance using the substrates complex index of refraction---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
						
						//---Now some wave clean up---//
						KillWaves/Z Nsub, DummyWave1, DummyWave2
					else
						//---Calculate the conplex conductance---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
					endif
					
					//---Sort waves into the appropriate folder---//
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Transmissions:$fname3
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Transmissions:$fname4
					
					//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
					Wave/C Dummy_Cond2 = $fname4 
					Dummy_Cond2 *= -sqrt(-1)
					
					//---Grabs the real and imaginary parts of the conductance---//
					Redimension/R $fname3
					Redimension/R $fname4
					
					//---Stores them in the appropriate folder---//
					Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Conductance1:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":Out_of_Phase:Conductance2:$fname4
					
					//---Some more wave clean up---//
					KillWaves/Z $fname2, $fname3, $fname4, freq
				
				endfor
			endif
endif

//---This is the case that we're running the Magnet---//
ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
	if (V_Value == 1)
	
		//---Define some variables---//
		variable NumFields, p;
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis":
			
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
			
				//---Creates folders for the in phase complex conductance of the thin film---//
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Conductance1
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Conductance2
				
				//---String for the in phase transmission---//
				fname = "Txx_" + Temps[i] + "_" + Fields[j] + "kG"
				
				//---String for the in phase conductance---//
				fname2 = ("Gxx_" + Temps[i] + "_" + Fields[j] + "kG")
				fname3 = ("G1xx_" + Temps[i] + "_" + Fields[j] + "kG")
				fname4 = ("G2xx_" + Temps[i] + "_" + Fields[j] + "kG")
				
				//---Now go to the transmissions folder to copy the waves to make conductance waves---//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:'Transmissions'
				Duplicate/O $Fname, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions:$fname2
				
				//---Makes a wave for the frequency--//
				delta = DimDelta($fname,0)
				Make/O/N=(numpnts($fname)), freq
				for(p = 0; p<numpnts(freq); p+=1)
					freq[p] = delta*p
				endfor
				SetScale/P x 0,(delta),"", freq
				
				//----dummy waves for calculations---//
				Wave/C Cond = $fname2
				Wave/C Trans = $fname
					
				//---Calculate the Complex Conductance---//
				ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
				If (V_Value == 1)
				
					//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_InField_Analysis":In_Phase:$Temps[i]:n
					Duplicate/O $("nxx_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions:$("nxx_" + Temps[i] + "_" + Fields[j] + "kG")
					
					//---Do the same thing with k, put it in the Sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_InField_Analysis":In_Phase:$Temps[i]:k
					Duplicate/O $("kxx_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions:$("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
					
					//---Now set the folder to the sample's transmissions's folder----//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:'Transmissions'
					
					//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
					Duplicate/O $fname, NSub
					
					//---Combine the n and k of the substrate into a single complex wave called NSub---//
					Wave/C Dummywave = Nsub
					Wave Dummywave1 = $("nxx_" + Temps[i] + "_" + Fields[j] + "kG")
					Wave Dummywave2 = $("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
					
					Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
				
					//---Calculate the conplex conductance using the substrates complex index of refraction---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
					
					//---Now some wave clean up---//
					KillWaves/Z Nsub, DummyWave1, DummyWave2
				else
					//---Calculate the conplex conductance---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
				endif
				
				//---Sort waves into the appropriate folder---//
				Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions:$fname3
				Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions:$fname4
				
				//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
				Wave/C Dummy_Cond2 = $fname4 
				Dummy_Cond2 *= -sqrt(-1)
				
				//---Grabs the real and imaginary parts of the conductance---//
				Redimension/R $fname3
				Redimension/R $fname4
				
				//---Stores them in the appropriate folder---//
				Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Conductance1:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Conductance2:$fname4
				
				//---Some more wave clean up---//
				KillWaves/Z $fname2, $fname3, $fname4, freq
			endfor
				
			//---This is the case where we're using the rotator---//
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			if (V_Value == 1)
				
				//---Creates folders for the out of phase complex conductance of the thin film---//
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Conductance1
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Conductance2
	
				for (j=0; j<NumFields; j+=1)
	
					//---String for the in phase transmission---//
					fname = "Tyx_" + Temps[i] + "_" + Fields[j] + "kG"
					
					//---String for the in phase conductance---//
					fname2 = ("Gyx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname3 = ("G1yx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname4 = ("G2yx_" + Temps[i] + "_" + Fields[j] + "kG")
					
					//---Now go to the transmissions folder to copy the waves to make conductance waves---//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:'Transmissions'
					Duplicate/O $Fname, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions:$fname2
					
					//---Makes a wave for the frequency--//
					delta = DimDelta($fname,0)
					Make/O/N=(numpnts($fname)), freq
					for(p = 0; p<numpnts(freq); p+=1)
						freq[p] = delta*p
					endfor
					SetScale/P x 0,(delta),"", freq
					
					//----dummy waves for calculations---//
					Wave/C Cond = $fname2
					Wave/C Trans = $fname
					
					//---Calculate the Complex Conductance---//
					ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
					If (V_Value == 1)
					
						//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n
						Duplicate/O $("nyx_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions:$("nyx_" + Temps[i] + "_" + Fields[j] + "kG")
						
						//---Do the same thing with k, put it in the Sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:k
						Duplicate/O $("kyx_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions:$("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
						
						//---Now set the folder to the sample's transmissions's folder----//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:'Transmissions'
						
						//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
						Duplicate/O $fname, NSub
						
						//---Combine the n and k of the substrate into a single complex wave called NSub---//
						Wave/C Dummywave = Nsub
						Wave Dummywave1 = $("nyx_" + Temps[i] + "_" + Fields[j] + "kG")
						Wave Dummywave2 = $("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
					
						//---Calculate the conplex conductance using the substrates complex index of refraction---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
						
						//---Now some wave clean up---//
						KillWaves/Z Nsub, DummyWave1, DummyWave2
					else
						//---Calculate the conplex conductance---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
					endif
					
					//---Sort waves into the appropriate folder---//
					Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions:$fname3
					Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions:$fname4
					
					//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
					Wave/C Dummy_Cond2 = $fname4 
					Dummy_Cond2 *= -sqrt(-1)
					
					//---Grabs the real and imaginary parts of the conductance---//
					Redimension/R $fname3
					Redimension/R $fname4
					
					//---Stores them in the appropriate folder---//
					Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Conductance1:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Conductance2:$fname4
					
					//---Some more wave clean up---//
					KillWaves/Z $fname2, $fname3, $fname4, freq
				
				endfor
			endif
		endfor
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the conductivity for a thin film for either the BLC or Magnet in the circular basis.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function CircBasisTFCond()

//---Folder where data is loaded from---//
SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//---String where you enter in the name of the sample, aperature and sub if chosen---//
SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
SVAR Sub_Name = root:Packages:NickLTHzExtras:Sub_Name

//---Substrate Index of Refraction, only used if the "use calculated substrate index of refraction" option is NOT chosen---//
NVAR ns = root:Packages:NickLTHzExtras:ns

//---DeltaL is difference in thicnkess between reference and sample substrate in microns----//
NVAR DeltaL = root:Packages:NickLTHzExtras:DeltaL

//---Dummy variables---//
String fname, fname2, fname3, fname4, fname5, fname6;
Variable i, j;
Variable Zo = 376.7, c = 3

//---Checks to see if we've already created the circular basis folder.
if (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:Circular_Basis") == 1)

	//---This is the case that we're running the BLC---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
	if (V_Value == 1)
	
		//---Sets the data folder to the Analysis folder just in case---//
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
		
		//---Recognizes the wave that we'll need while running this function---//
		Wave/T Temps
		Wave num_samp_scans
		Variable NumTemps = numpnts(Temps)
		
			//---Creates folders for the complex conductance of the thin film---//
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Conductance1
			NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Conductance2
				
			for(i = 0; i < NumTemps; i+=1)
				
				//--Strings for the in phase transmissions--//
				fname = "Tr_" + Sample_Name + "_" + Temps[i]
				
				//--Strings for the in phase Conductance--//
				fname2 = "Gr_" + Sample_Name + "_"+ Temps[i]
				fname3 =  "G1r_" + Sample_Name + "_"+ Temps[i]
				fname4 =  "G2r_" + Sample_Name + "_"+ Temps[i]
				
				//---Now go to the transmissions folder to copy the waves to make conductance waves---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:'Transmissions'
				Duplicate/O $fname, $fname2
				Duplicate/O $fname, $fname3
				Duplicate/O $fname, $fname4
				
				Redimension/R $fname3, $fname4
				
				//---Makes a wave for the frequency--//
				Variable delta = DimDelta($fname,0)
				Make/O/N=(numpnts($fname)), freq
				for(j = 0; j<numpnts(freq); j+=1)
					freq[j] = delta*j
				endfor
				SetScale/P x 0,(delta),"", freq
				
				//---dummy waves for calculations---//
				Wave/C Trans = $fname
				Wave/C Cond = $fname2
				Wave Cond1 = $fname3
				Wave Cond2 = $fname4
				
				//---This is the option that we're using an already calculated index of refraction of the substrate via this analysis palette---//
				//---The analysis folder of the substrate has to exist prior to trying to analyze the thin films.  So first run the entire code on the sub as a single crystal---//
				//---Then run the code on the thin film---//
				ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
				If (V_Value == 1)
					
					//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n
					Duplicate/O $("nr_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Transmissions:$("nr_" + Sub_Name + "_" + Temps[i])
					
					//---Do the same thing with k, put it in the Sample's transmissions folder---//
					SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:k
					Duplicate/O $("kr_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Transmissions:$("kr_" + Sub_Name + "_" + Temps[i])
					
					//---Now set the folder to the sample's transmissions's folder----//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:'Transmissions'
					
					//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
					Duplicate/O $fname, NSub
					
					//---Combine the n and k of the substrate into a single complex wave called NSub---//
					Wave/C Dummywave = Nsub
					Wave Dummywave1 = $("nr_" + Sub_Name + "_" + Temps[i])
					Wave Dummywave2 = $("kr_" + Sub_Name + "_" + Temps[i])
					
					Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
				
					//---Calculate the conplex conductance using the substrates complex index of refraction---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
					
					//---Now some wave clean up---//
					KillWaves/Z Nsub, DummyWave1, DummyWave2
				else
					//---Calculate the conplex conductance---//
					cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
				endif
				
				//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
				Cond1 = Real(Cond)
				Cond2 = Imag(Cond)

				//---Stores them in the appropriate folder---//
				Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Conductance1:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Conductance2:$fname4
				
				//---Some more wave clean up---//
				KillWaves/Z $fname2, $fname3, $fname4, freq
			
			endfor
			
			//---Now repeat the entire process for the out of phase components---//
			//---This is the case that we're using the rotator
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			if (V_Value == 1)
				
				//---Creates folders for the complex conductance of the thin film---//
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Conductance1
				NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Conductance2
					
				for(i = 0; i < numpnts(Temps); i+=1)
					
					//--Strings for the in phase transmissions--//
					fname = "Tl_" + Sample_Name + "_" + Temps[i]
					
					//--Strings for the in phase Conductance--//
					fname2 = "Gl_" + Sample_Name + "_"+ Temps[i]
					fname3 =  "G1l_" + Sample_Name + "_"+ Temps[i]
					fname4 =  "G2l_" + Sample_Name + "_"+ Temps[i]
					
					//---Now go to the transmissions folder to copy the waves to make conductance waves---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:'Transmissions'
					Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$fname2
					
					//---Makes a wave for the frequency--//
					delta = DimDelta($fname,0)
					Make/O/N=(numpnts($fname)), freq
					for(j = 0; j<numpnts(freq); j+=1)
						freq[j] = delta*j
					endfor
					SetScale/P x 0,(delta),"", freq
					
					//----dummy waves for calculations---//
					Wave/C Cond = $fname2
					Wave/C Trans = $fname
					
					//---Calculate the Complex Conductance---//
					ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
					If (V_Value == 1)
					
						//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n
						Duplicate/O $("nl_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$("nl_" + Sub_Name + "_" + Temps[i])
						
						//---Do the same thing with k, put it in the Sample's transmissions folder---//
						SetDataFolder  root:$Sub_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:k
						Duplicate/O $("kl_" + Sub_Name + "_" + Temps[i]),  root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$("kl_" + Sub_Name + "_" + Temps[i])
						
						//---Now set the folder to the sample's transmissions's folder----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:'Transmissions'
						
						//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
						Duplicate/O $fname, NSub
						
						//---Combine the n and k of the substrate into a single complex wave called NSub---//
						Wave/C Dummywave = Nsub
						Wave Dummywave1 = $("nl_" + Sub_Name + "_" + Temps[i])
						Wave Dummywave2 = $("kl_" + Sub_Name + "_" + Temps[i])
						
						Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
					
						//---Calculate the conplex conductance using the substrates complex index of refraction---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
						
						//---Now some wave clean up---//
						KillWaves/Z Nsub, DummyWave1, DummyWave2
					else
						//---Calculate the conplex conductance---//
						cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
					endif
					
					//---Sort waves into the appropriate folder---//
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$fname3
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions:$fname4
					
					//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
					Wave/C Dummy_Cond2 = $fname4 
					Dummy_Cond2 *= -sqrt(-1)
					
					//---Grabs the real and imaginary parts of the conductance---//
					Redimension/R $fname3
					Redimension/R $fname4
					
					//---Stores them in the appropriate folder---//
					Duplicate/O $fname3, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Conductance1:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Conductance2:$fname4
					
					//---Some more wave clean up---//
					KillWaves/Z $fname2, $fname3, $fname4, freq
				
				endfor
			endif
	endif
endif

//---Checks to see if we've already created the circular basis folder.
if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:Circular_Basis") == 1)

		//---This is the case that we're running the Magnet---//
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
			if (V_Value == 1)
			
				//---Define some variables---//
				variable NumFields, p;
			
				SetDataFolder root:$Sample_Name + "_InField_Analysis":
				Wave/T Temps
				NumTemps = NumPnts(Temps)
				
				for (i=0; i<NumTemps; i+=1)
				
					//---Set the data folder to the main folder to find "Fields_Temps" waves
					SetDataFolder root:$Sample_Name + "_InField_Analysis":
					
					fname = "Fields_" + Temps[i]
					Wave/T Fields = $fname
					NumFields = NumPnts($fname)
				
					//--- j runs over the number of Fields we have---//
					for (j=0; j<NumFields; j+=1)
					
						//---Creates folders for the in phase complex conductance of the thin film---//
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Conductance1
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Conductance2
						
						//---String for the in phase transmission---//
						fname = "Tr_" + Temps[i] + "_" + Fields[j] + "kG"
						
						//---String for the in phase conductance---//
						fname2 = ("Gr_" + Temps[i] + "_" + Fields[j] + "kG")
						fname3 = ("G1r_" + Temps[i] + "_" + Fields[j] + "kG")
						fname4 = ("G2r_" + Temps[i] + "_" + Fields[j] + "kG")
						
						//---Now go to the transmissions folder to copy the waves to make conductance waves---//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:'Transmissions'
						Duplicate/O $Fname, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$fname2
						
						//---Makes a wave for the frequency--//
						delta = DimDelta($fname,0)
						Make/O/N=(numpnts($fname)), freq
						for(p = 0; p<numpnts(freq); p+=1)
							freq[p] = delta*p
						endfor
						SetScale/P x 0,(delta),"", freq
						
						//----dummy waves for calculations---//
						Wave/C Cond = $fname2
						Wave/C Trans = $fname
							
						//---Calculate the Complex Conductance---//
						ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
						If (V_Value == 1)
						
							//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
							SetDataFolder  root:$Sub_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n
							Duplicate/O $("nr_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$("nr_" + Temps[i] + "_" + Fields[j] + "kG")
							
							//---Do the same thing with k, put it in the Sample's transmissions folder---//
							SetDataFolder  root:$Sub_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:k
							Duplicate/O $("kr_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$("kr_" + Temps[i] + "_" + Fields[j] + "kG")
							
							//---Now set the folder to the sample's transmissions's folder----//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:'Transmissions'
							
							//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
							Duplicate/O $fname, NSub
							
							//---Combine the n and k of the substrate into a single complex wave called NSub---//
							Wave/C Dummywave = Nsub
							Wave Dummywave1 = $("nr_" + Temps[i] + "_" + Fields[j] + "kG")
							Wave Dummywave2 = $("kr_" + Temps[i] + "_" + Fields[j] + "kG")
							
							Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
						
							//---Calculate the conplex conductance using the substrates complex index of refraction---//
							cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
							
							//---Now some wave clean up---//
							KillWaves/Z Nsub, DummyWave1, DummyWave2
						else
							//---Calculate the conplex conductance---//
							cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
						endif
						
						//---Sort waves into the appropriate folder---//
						Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$fname3
						Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions:$fname4
						
						//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
						Wave/C Dummy_Cond2 = $fname4 
						Dummy_Cond2 *= -sqrt(-1)
						
						//---Grabs the real and imaginary parts of the conductance---//
						Redimension/R $fname3
						Redimension/R $fname4
						
						//---Stores them in the appropriate folder---//
						Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Conductance1:$fname3
						Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Conductance2:$fname4
						
						//---Some more wave clean up---//
						KillWaves/Z $fname2, $fname3, $fname4, freq
					endfor
						
					//---This is the case where we're using the rotator---//
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)
						
						//---Creates folders for the out of phase complex conductance of the thin film---//
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Conductance1
						NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Conductance2
			
						for (j=0; j<NumFields; j+=1)
			
							//---String for the in phase transmission---//
							fname = "Tl_" + Temps[i] + "_" + Fields[j] + "kG"
							
							//---String for the in phase conductance---//
							fname2 = ("Gl_" + Temps[i] + "_" + Fields[j] + "kG")
							fname3 = ("G1l_" + Temps[i] + "_" + Fields[j] + "kG")
							fname4 = ("G2l_" + Temps[i] + "_" + Fields[j] + "kG")
							
							//---Now go to the transmissions folder to copy the waves to make conductance waves---//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:'Transmissions'
							Duplicate/O $Fname, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$fname2
							
							//---Makes a wave for the frequency--//
							delta = DimDelta($fname,0)
							Make/O/N=(numpnts($fname)), freq
							for(p = 0; p<numpnts(freq); p+=1)
								freq[p] = delta*p
							endfor
							SetScale/P x 0,(delta),"", freq
							
							//----dummy waves for calculations---//
							Wave/C Cond = $fname2
							Wave/C Trans = $fname
							
							//---Calculate the Complex Conductance---//
							ControlInfo/W=NickL_THzAnalysis_Palette Subn_CheckBox
							If (V_Value == 1)
							
								//---First go to the Sub folder and grab n and duplicate it into the sample's transmissions folder---//
								SetDataFolder  root:$Sub_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n
								Duplicate/O $("nl_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$("nl_" + Temps[i] + "_" + Fields[j] + "kG")
								
								//---Do the same thing with k, put it in the Sample's transmissions folder---//
								SetDataFolder  root:$Sub_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:k
								Duplicate/O $("kl_" + Temps[i] + "_" + Fields[j] + "kG"),  root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$("kl_" + Temps[i] + "_" + Fields[j] + "kG")
								
								//---Now set the folder to the sample's transmissions's folder----//
								SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:'Transmissions'
								
								//---Duplicate the transmission so we have a wave for the complex substrate index of refraction, call it Nsub---//
								Duplicate/O $fname, NSub
								
								//---Combine the n and k of the substrate into a single complex wave called NSub---//
								Wave/C Dummywave = Nsub
								Wave Dummywave1 = $("nl_" + Temps[i] + "_" + Fields[j] + "kG")
								Wave Dummywave2 = $("kl_" + Temps[i] + "_" + Fields[j] + "kG")
								
								Dummywave = Dummywave1 + sqrt(-1)*Dummywave2
							
								//---Calculate the conplex conductance using the substrates complex index of refraction---//
								cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(NSub-1)*0.01*(DeltaL))-1)*(NSub+1)/Zo
								
								//---Now some wave clean up---//
								KillWaves/Z Nsub, DummyWave1, DummyWave2
							else
								//---Calculate the conplex conductance---//
								cond=(1/Trans*exp(-1*sqrt(-1)*2*pi*freq/c*(ns-1)*0.01*(DeltaL))-1)*(ns+1)/Zo
							endif
							
							//---Sort waves into the appropriate folder---//
							Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$fname3
							Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions:$fname4
							
							//---We want to split complex conductance into the real and imaginary parts.  We make a dummy wave to store the imaginary part---//
							Wave/C Dummy_Cond2 = $fname4 
							Dummy_Cond2 *= -sqrt(-1)
							
							//---Grabs the real and imaginary parts of the conductance---//
							Redimension/R $fname3
							Redimension/R $fname4
							
							//---Stores them in the appropriate folder---//
							Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Conductance1:$fname3
							Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Conductance2:$fname4
							
							//---Some more wave clean up---//
							KillWaves/Z $fname2, $fname3, $fname4, freq
						
						endfor
					endif
				endfor
			endif
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the Susceptibility for a single crystal for either the BLC or Magnet, with or without the rotator.  Grabs an index of refraction at some higher temperature above the magnetism in the sample which
//  the user specifies for the calculation.  This temperature is call Ref_Temp.  Note there is no function for calculating the susceptibility of a thin film.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
Function Susceptibility()

//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name

//-----Thickness of the sample in mm's------//
	NVAR dsc = root:Packages:NickLTHzExtras:dsc

//-----PopVal sets the type of sample being measured so we know which version of NR code to run------//
	NVAR PopVal = root:Packages:NickLTHzExtras:PopVal

//-----Substrate Index of Refraction------//
	NVAR ns = root:Packages:NickLTHzExtras:ns

//-----Reference Temperature for Susceptibility Calculation------//
	SVAR RefTemp = root:Packages:NickLTHzExtras:RefTemp
	
//------Declare other variables and dummy waves------//
	variable delta, i, j, c = 299792458, mu = 1//(1.25663706*10^-6);
	Variable num_iter = 10  //--The number of Netwon-Rhapson iterations
	
//-----Dummy variables-----//
	String fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9
	
	//---This is the case that we're running the BLC---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		if (V_Value == 1)
		
			//---Sets the data folder to the Analysis folder just in case---//
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
			
			//---recognizes the wave that we'll need while running this function---//
			Wave/T Temps
			Wave num_samp_scans
			 Variable NumTemps = numpnts(Temps)
			
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:k
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Phase
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Chi1
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Chi2
					
					for(i = 0; i < NumTemps; i+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Transmissions
						index_calc("Txx_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("Txx_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $("n"), root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n:$("nxx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("n")
						
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n:$("kxx_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":In_Phase:k:$("kxx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("k")
						
						Duplicate/O $("phase"), root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Phase:$("Thetaxx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("phase")
						
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					
					for(i = 0; i < numpnts(Temps); i+=1)
						
						//---Now go to the Mag Transmissions folder so we can calculate the complex susceptibility.  n, k, phase, and MagT should all be there already---//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:n
						
						//---Name some strings for waves we'll need later---//
						fname1 = ("nxx_" + Sample_Name + "_" + Temps[i])
						fname2 = ("kxx_" + Sample_Name + "_" + Temps[i])
						
						fname3 = ("nxx_" + Sample_Name + "_" + RefTemp)
						fname4 = ("kxx_" + Sample_Name + "_" + RefTemp)
						
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname6 = ("nComplex_" + Sample_Name + "_" + RefTemp)
						
						fname7 = ("Chi2xx_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1xx_" + Sample_Name + "_" + Temps[i])
						
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						//---We'll need a wave of omega, the angular frequency---//
						Make/C/O/N=(numpnts($fname1)), $fname5
						Make/C/O/N=(numpnts($fname1)), $fname6
						Make/O/N=(numpnts($fname1)), $fname7
						Make/O/N=(numpnts($fname1)), $fname8
						Make/C/O/N=(numpnts($fname1)), $fname9
						
						//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
						Wave F1 = $fname1
						Wave F2 = $fname2
						Wave F3 = $fname3
						Wave F4 = $fname4
						Wave/C F5 = $fname5
						Wave/C F6 = $fname6
						Wave F7 = $fname7
						Wave F8 = $fname8
						Wave/C F9 = $fname9
					
						delta = DimDelta($fname1,0)
						
						SetScale/P x 0,(delta),"", $fname5
						SetScale/P x 0,(delta),"", $fname6
						SetScale/P x 0,(delta),"", $fname7
						SetScale/P x 0,(delta),"", $fname8
						SetScale/P x 0,(delta),"", $fname9
						
						F6 = F3 //+ sqrt(-1)*F4
						F5 = F1 + sqrt(-1)*F2
						
						F9 = (F5/F6)^2 - 1
						
						F8 = real(F9)
						F7 = imag(F9)
						
					endfor
					
					//---A loop which will do some wave clean up for us---//
					for(i = 0; i < numpnts(Temps); i+=1)
						
						fname2 = ("kxx_" + Sample_Name + "_" + Temps[i])
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname7 = ("Chi2xx_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1xx_" + Sample_Name + "_" + Temps[i])
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						Duplicate/O $fname7, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
					endfor
					
			//---This is the case that we're using the rotator---//
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
				if (V_Value == 1)
					
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:k
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Phase
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Chi1
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Chi2
					
					for(i = 0; i < numpnts(Temps); i+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Transmissions
						index_calc("Tyx_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("Tyx_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $("n"), root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n:$("nyx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("n")
						
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n:$("kyx_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:k:$("kyx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("k")
						
						Duplicate/O $("phase"), root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Phase:$("Thetayx_" + Sample_Name + "_" + Temps[i])
						KillWaves $("phase")
						
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					
					for(i = 0; i < numpnts(Temps); i+=1)
						
						//---Now go to the Mag Transmissions folder so we can calculate the complex susceptibility.  n, k, phase, and MagT should all be there already---//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:n
						
						//---Name some strings for waves we'll need later---//
						fname1 = ("nyx_" + Sample_Name + "_" + Temps[i])
						fname2 = ("kyx_" + Sample_Name + "_" + Temps[i])
						
						fname3 = ("nyx_" + Sample_Name + "_" + RefTemp)
						fname4 = ("kyx_" + Sample_Name + "_" + RefTemp)
						
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname6 = ("nComplex_" + Sample_Name + "_" + RefTemp)
						
						fname7 = ("Chi2yx_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1yx_" + Sample_Name + "_" + Temps[i])
						
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						//---We'll need a wave of omega, the angular frequency---//
						Make/C/O/N=(numpnts($fname1)), $fname5
						Make/C/O/N=(numpnts($fname1)), $fname6
						Make/O/N=(numpnts($fname1)), $fname7
						Make/O/N=(numpnts($fname1)), $fname8
						Make/C/O/N=(numpnts($fname1)), $fname9
						
						//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
						Wave F1 = $fname1
						Wave F2 = $fname2
						Wave F3 = $fname3
						Wave F4 = $fname4
						Wave/C F5 = $fname5
						Wave/C F6 = $fname6
						Wave F7 = $fname7
						Wave F8 = $fname8
						Wave/C F9 = $fname9
					
						delta = DimDelta($fname1,0)
						
						SetScale/P x 0,(delta),"", $fname5
						SetScale/P x 0,(delta),"", $fname6
						SetScale/P x 0,(delta),"", $fname7
						SetScale/P x 0,(delta),"", $fname8
						SetScale/P x 0,(delta),"", $fname9
						
						F6 = F3 //+ sqrt(-1)*F4
						F5 = F1 + sqrt(-1)*F2
						
						F9 = (F5/F6)^2 - 1
						
						F8 = real(F9)
						F7 = imag(F9)
						
					endfor
					
					//---A loop which will do some wave clean up for us---//
					for(i = 0; i < numpnts(Temps); i+=1)
						
						fname2 = ("kyx_" + Sample_Name + "_" + Temps[i])
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname7 = ("Chi2yx_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1yx_" + Sample_Name + "_" + Temps[i])
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						Duplicate/O $fname7, root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
					endfor
			endif
		endif
		
	//---This is the case that we're running the Magnet---//
	ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
	if (V_Value == 1)
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis"
		
			//---Creates folders for the complex conductivity, complex n, and phase---//
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:k
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Phase
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi1
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi2
			
			string fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			Variable NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
		
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions
				index_calc("Txx_" + Temps[i] + "_" + Fields[j] + "kG")
				Wave transmission, phase, freq
				
				//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
				CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
				Wave W_coef
				phase -= W_coef[0]
				
				//-----Run the newton-raphson code to find the index-----//
				delta = DimDelta($("Txx_" + Temps[i] + "_" + Fields[j] + "kG"),0)
				
				//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
				//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
				if (PopVal==1)
					//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
					TDNR(transmission, phase, freq, dsc, num_iter)
				elseif (PopVal==2)
					//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
					TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
				endif
				
				//---Declare n and k and then set their scale---//
				Wave n,k
				SetScale/P x 0,(delta),"", n,k
				
				//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
				Duplicate/O $("n"), root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n:$("nxx_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("n")
				
				Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n:$("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:k:$("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("k")
				
				Duplicate/O $("phase"), root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Phase:$("Thetaxx_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("phase")
					
			endfor
				
					//-----A bit of wave clean up-----//
					Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
		endfor
		
		for (i=0; i<NumTemps; i+=1)	
			SetDataFolder root:$Sample_Name + "_InField_Analysis"
					
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
			
			for (j=0; j<NumFields; j+=1)
				
				//---Name some strings for waves we'll need later---//
				fname1 = ("nxx_" + Temps[i] + "_" + Fields[j] + "kG")
				fname2 = ("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
				fname3 = ("nxx_" + RefTemp + "_" + Fields[j] + "kG")
				fname4 = ("kxx_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
				fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname7 = ("Chi2xx_" + Temps[i] + "_" + Fields[j] + "kG")
				fname8 = ("Chi1xx_" +  Temps[i] + "_" + Fields[j] + "kG")
				
				fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
				
				SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$(RefTemp):n
				Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n:$fname4
				
				SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:n
				
				//---Make some waves we'll need for the calculation---//
				Make/C/O/N=(numpnts($fname1)), $fname5
				Make/C/O/N=(numpnts($fname1)), $fname6
				Make/O/N=(numpnts($fname1)), $fname7
				Make/O/N=(numpnts($fname1)), $fname8
				Make/C/O/N=(numpnts($fname1)), $fname9
				
				//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
				Wave F1 = $fname1
				Wave F2 = $fname2
				Wave F3 = $fname3
				Wave F4 = $fname4
				Wave/C F5 = $fname5
				Wave/C F6 = $fname6
				Wave F7 = $fname7
				Wave F8 = $fname8
				Wave/C F9 = $fname9
			
				delta = DimDelta($fname1,0)
				
				SetScale/P x 0,(delta),"", $fname5
				SetScale/P x 0,(delta),"", $fname6
				SetScale/P x 0,(delta),"", $fname7
				SetScale/P x 0,(delta),"", $fname8
				SetScale/P x 0,(delta),"", $fname9
				
				F6 = F3 //+ sqrt(-1)*F4
				F5 = F1 + sqrt(-1)*F2
				
				F9 = (F5/F6)^2 - 1
				
				F8 = real(F9)
				F7 = imag(F9)
		
			endfor
				
				//---A loop which will do some wave clean up for us---//
				
				if (cmpstr(Temps[i], RefTemp) == 0)
					for (j=0; j<NumFields; j+=1)
					
						fname2 = ("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname7 = ("Chi2xx_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1xx_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
						
					endfor
					
				else
					for (j=0; j<NumFields; j+=1)
						fname2 = ("kxx_" + Temps[i] + "_" + Fields[j] + "kG")
						fname3 = ("nxx_" + RefTemp + "_" + Fields[j] + "kG")
						fname4 = ("kxx_" + RefTemp + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
						fname7 = ("Chi2xx_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1xx_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname3, $fname4, $fname5, $fname6, $fname7, $fname8, $fname9
					endfor
				endif
			endfor
		
		//---This is the case that we're using the rotator---//
		ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
		if (V_Value == 1)
		
			for (i=0; i<NumTemps; i+=1)
			
				//---Set the data folder to the main folder to find "Fields_Temps" waves
				SetDataFolder root:$Sample_Name + "_InField_Analysis"
			
				//---Creates folders for the complex conductivity, complex n, and phase---//
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:k
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Phase
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi1
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi2
				
				fname = "Fields_" + Temps[i]
				Wave/T Fields = $fname
				NumFields = NumPnts($fname)
			
				//--- j runs over the number of Fields we have---//
				for (j=0; j<NumFields; j+=1)
			
					//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Transmissions
					index_calc("Tyx_" + Temps[i] + "_" + Fields[j] + "kG")
					Wave transmission, phase, freq
					
					//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
					CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
					Wave W_coef
					phase -= W_coef[0]
					
					//-----Run the newton-raphson code to find the index-----//
					delta = DimDelta($("Tyx_" + Temps[i] + "_" + Fields[j] + "kG"),0)
					
					//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
					//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
					if (PopVal==1)
						//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
						TDNR(transmission, phase, freq, dsc, num_iter)
					elseif (PopVal==2)
						//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
						TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
					endif
					
					//---Declare n and k and then set their scale---//
					Wave n,k
					SetScale/P x 0,(delta),"", n,k
					
					//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
					Duplicate/O $("n"), root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n:$("nyx_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("n")
					
					Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n:$("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:k:$("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("k")
					
					Duplicate/O $("phase"), root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Phase:$("Thetayx_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("phase")
						
				endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			endfor
			
			for (i=0; i<NumTemps; i+=1)	
				SetDataFolder root:$Sample_Name + "_InField_Analysis"
						
				fname = "Fields_" + Temps[i]
				Wave/T Fields = $fname
				NumFields = NumPnts($fname)
				
				for (j=0; j<NumFields; j+=1)
					
					//---Name some strings for waves we'll need later---//
					fname1 = ("nyx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname2 = ("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname3 = ("nyx_" + RefTemp + "_" + Fields[j] + "kG")
					fname4 = ("kyx_" + RefTemp + "_" + Fields[j] + "kG")
					
					fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
					fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
					
					fname7 = ("Chi2yx_" + Temps[i] + "_" + Fields[j] + "kG")
					fname8 = ("Chi1yx_" +  Temps[i] + "_" + Fields[j] + "kG")
					
					fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
					
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$(RefTemp):n
					Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n:$fname4
					
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:n
					
					//---Make some waves we'll need for the calculation---//
					Make/C/O/N=(numpnts($fname1)), $fname5
					Make/C/O/N=(numpnts($fname1)), $fname6
					Make/O/N=(numpnts($fname1)), $fname7
					Make/O/N=(numpnts($fname1)), $fname8
					Make/C/O/N=(numpnts($fname1)), $fname9
					
					//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
					Wave F1 = $fname1
					Wave F2 = $fname2
					Wave F3 = $fname3
					Wave F4 = $fname4
					Wave/C F5 = $fname5
					Wave/C F6 = $fname6
					Wave F7 = $fname7
					Wave F8 = $fname8
					Wave/C F9 = $fname9
				
					delta = DimDelta($fname1,0)
					
					SetScale/P x 0,(delta),"", $fname5
					SetScale/P x 0,(delta),"", $fname6
					SetScale/P x 0,(delta),"", $fname7
					SetScale/P x 0,(delta),"", $fname8
					SetScale/P x 0,(delta),"", $fname9
					
					F6 = F3 //+ sqrt(-1)*F4
					F5 = F1 + sqrt(-1)*F2
					
					F9 = (F5/F6)^2 - 1
					
					F8 = real(F9)
					F7 = imag(F9)
			
				endfor
					
					//---A loop which will do some wave clean up for us---//
					
					if (cmpstr(Temps[i], RefTemp) == 0)
						for (j=0; j<NumFields; j+=1)
						
							fname2 = ("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
							fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
							fname7 = ("Chi2yx_" + Temps[i] + "_" + Fields[j] + "kG")
							fname8 = ("Chi1yx_" +  Temps[i] + "_" + Fields[j] + "kG")
							fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
							
							Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi2:$fname7
							Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi1:$fname8
							
							Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
							
						endfor
						
					else
						for (j=0; j<NumFields; j+=1)
							fname2 = ("kyx_" + Temps[i] + "_" + Fields[j] + "kG")
							fname3 = ("nyx_" + RefTemp + "_" + Fields[j] + "kG")
							fname4 = ("kyx_" + RefTemp + "_" + Fields[j] + "kG")
							fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
							fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
							fname7 = ("Chi2yx_" + Temps[i] + "_" + Fields[j] + "kG")
							fname8 = ("Chi1yx_" +  Temps[i] + "_" + Fields[j] + "kG")
							fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
							
							Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi2:$fname7
							Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Out_of_Phase:$Temps[i]:Chi1:$fname8
							
							Killwaves/Z $fname2, $fname3, $fname4, $fname5, $fname6, $fname7, $fname8, $fname9
						endfor
					endif
			endfor
		endif
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the Susceptibility for a single crystal for either the BLC or Magnet in the Circular Basis.  Grabs an index of refraction at some higher temperature above the magnetism in the sample which
//  the user specifies for the calculation.  This temperature is call Ref_Temp.  Note there is no function for calculating the susceptibility of a thin film.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
Function CircBasisSusceptibility()

//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name

//-----Thickness of the sample in mm's------//
	NVAR dsc = root:Packages:NickLTHzExtras:dsc

//-----PopVal sets the type of sample being measured so we know which version of NR code to run------//
	NVAR PopVal = root:Packages:NickLTHzExtras:PopVal

//-----Substrate Index of Refraction------//
	NVAR ns = root:Packages:NickLTHzExtras:ns

//-----Reference Temperature for Susceptibility Calculation------//
	SVAR RefTemp = root:Packages:NickLTHzExtras:RefTemp
	
//------Declare other variables and dummy waves------//
	variable delta, i, j, c = 299792458, mu = 1//(1.25663706*10^-6);
	Variable num_iter = 10  //--The number of Netwon-Rhapson iterations
	
//-----Dummy variables-----//
	String fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9
	
//---Checks to see if we've already created the circular basis folder.
if (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:Circular_Basis") == 1)
	
	//---This is the case that we're running the BLC---//
	ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		if (V_Value == 1)
		
			//---Sets the data folder to the Analysis folder just in case---//
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
			
			//---recognizes the wave that we'll need while running this function---//
			Wave/T Temps
			Wave num_samp_scans
			 Variable NumTemps = numpnts(Temps)
			
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:k
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Phase
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Chi1
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Chi2
					
					for(i = 0; i < NumTemps; i+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Transmissions
						index_calc("Tr_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("Tr_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $("n"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n:$("nr_" + Sample_Name + "_" + Temps[i])
						KillWaves $("n")
						
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n:$("kr_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:k:$("kr_" + Sample_Name + "_" + Temps[i])
						KillWaves $("k")
						
						Duplicate/O $("phase"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Phase:$("Thetar_" + Sample_Name + "_" + Temps[i])
						KillWaves $("phase")
						
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					
					for(i = 0; i < numpnts(Temps); i+=1)
						
						//---Now go to the Mag Transmissions folder so we can calculate the complex susceptibility.  n, k, phase, and MagT should all be there already---//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:n
						
						//---Name some strings for waves we'll need later---//
						fname1 = ("nr_" + Sample_Name + "_" + Temps[i])
						fname2 = ("kr_" + Sample_Name + "_" + Temps[i])
						
						fname3 = ("nr_" + Sample_Name + "_" + RefTemp)
						fname4 = ("kr_" + Sample_Name + "_" + RefTemp)
						
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname6 = ("nComplex_" + Sample_Name + "_" + RefTemp)
						
						fname7 = ("Chi2r_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1r_" + Sample_Name + "_" + Temps[i])
						
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						//---We'll need a wave of omega, the angular frequency---//
						Make/C/O/N=(numpnts($fname1)), $fname5
						Make/C/O/N=(numpnts($fname1)), $fname6
						Make/O/N=(numpnts($fname1)), $fname7
						Make/O/N=(numpnts($fname1)), $fname8
						Make/C/O/N=(numpnts($fname1)), $fname9
						
						//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
						Wave F1 = $fname1
						Wave F2 = $fname2
						Wave F3 = $fname3
						Wave F4 = $fname4
						Wave/C F5 = $fname5
						Wave/C F6 = $fname6
						Wave F7 = $fname7
						Wave F8 = $fname8
						Wave/C F9 = $fname9
					
						delta = DimDelta($fname1,0)
						
						SetScale/P x 0,(delta),"", $fname5
						SetScale/P x 0,(delta),"", $fname6
						SetScale/P x 0,(delta),"", $fname7
						SetScale/P x 0,(delta),"", $fname8
						SetScale/P x 0,(delta),"", $fname9
						
						F6 = F3 //+ sqrt(-1)*F4
						F5 = F1 + sqrt(-1)*F2
						
						F9 = (F5/F6)^2 - 1
						
						F8 = real(F9)
						F7 = imag(F9)
						
					endfor
					
					//---A loop which will do some wave clean up for us---//
					for(i = 0; i < numpnts(Temps); i+=1)
						
						fname2 = ("kr_" + Sample_Name + "_" + Temps[i])
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname7 = ("Chi2r_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1r_" + Sample_Name + "_" + Temps[i])
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						Duplicate/O $fname7, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:RightHand:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
					endfor
					
			//---This is the case that we're using the rotator---//
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
				if (V_Value == 1)
					
					//---Creates folders for the complex conductivity, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:k
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Phase
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Chi1
					NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Chi2
					
					for(i = 0; i < numpnts(Temps); i+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Transmissions
						index_calc("Tl_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("Tl_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
						//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
						if (PopVal==1)
							//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
							TDNR(transmission, phase, freq, dsc, num_iter)
						elseif (PopVal==2)
							//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
							TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
						endif
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $("n"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n:$("nl_" + Sample_Name + "_" + Temps[i])
						KillWaves $("n")
						
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n:$("kl_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $("k"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:k:$("kl_" + Sample_Name + "_" + Temps[i])
						KillWaves $("k")
						
						Duplicate/O $("phase"), root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Phase:$("Thetal_" + Sample_Name + "_" + Temps[i])
						KillWaves $("phase")
						
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					
					for(i = 0; i < numpnts(Temps); i+=1)
						
						//---Now go to the Mag Transmissions folder so we can calculate the complex susceptibility.  n, k, phase, and MagT should all be there already---//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:n
						
						//---Name some strings for waves we'll need later---//
						fname1 = ("nl_" + Sample_Name + "_" + Temps[i])
						fname2 = ("kl_" + Sample_Name + "_" + Temps[i])
						
						fname3 = ("nl_" + Sample_Name + "_" + RefTemp)
						fname4 = ("kl_" + Sample_Name + "_" + RefTemp)
						
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname6 = ("nComplex_" + Sample_Name + "_" + RefTemp)
						
						fname7 = ("Chi2l_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1l_" + Sample_Name + "_" + Temps[i])
						
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						//---We'll need a wave of omega, the angular frequency---//
						Make/C/O/N=(numpnts($fname1)), $fname5
						Make/C/O/N=(numpnts($fname1)), $fname6
						Make/O/N=(numpnts($fname1)), $fname7
						Make/O/N=(numpnts($fname1)), $fname8
						Make/C/O/N=(numpnts($fname1)), $fname9
						
						//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
						Wave F1 = $fname1
						Wave F2 = $fname2
						Wave F3 = $fname3
						Wave F4 = $fname4
						Wave/C F5 = $fname5
						Wave/C F6 = $fname6
						Wave F7 = $fname7
						Wave F8 = $fname8
						Wave/C F9 = $fname9
					
						delta = DimDelta($fname1,0)
						
						SetScale/P x 0,(delta),"", $fname5
						SetScale/P x 0,(delta),"", $fname6
						SetScale/P x 0,(delta),"", $fname7
						SetScale/P x 0,(delta),"", $fname8
						SetScale/P x 0,(delta),"", $fname9
						
						F6 = F3 //+ sqrt(-1)*F4
						F5 = F1 + sqrt(-1)*F2
						
						F9 = (F5/F6)^2 - 1
						
						F8 = real(F9)
						F7 = imag(F9)
						
					endfor
					
					//---A loop which will do some wave clean up for us---//
					for(i = 0; i < numpnts(Temps); i+=1)
						
						fname2 = ("kl_" + Sample_Name + "_" + Temps[i])
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname7 = ("Chi2l_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1l_" + Sample_Name + "_" + Temps[i])
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						Duplicate/O $fname7, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_ZeroField_Analysis":Circular_Basis:LeftHand:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
					endfor
			endif
		endif
endif

//---Checks to see if we've already created the circular basis folder.
if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:Circular_Basis") == 1)

	//---This is the case that we're running the Magnet---//
	ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
	if (V_Value == 1)
	
	//---This is the case that we're using the rotator---//
	ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
	if (V_Value == 1)
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_InField_Analysis"
		
			//---Creates folders for the complex conductivity, complex n, and phase---//
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:k
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Phase
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi1
			NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi2
			
			string fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			Variable NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
		
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
				SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Transmissions
				index_calc("Tr_" + Temps[i] + "_" + Fields[j] + "kG")
				Wave transmission, phase, freq
				
				//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
				CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
				Wave W_coef
				phase -= W_coef[0]
				
				//-----Run the newton-raphson code to find the index-----//
				delta = DimDelta($("Tr_" + Temps[i] + "_" + Fields[j] + "kG"),0)
				
				//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
				//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
				if (PopVal==1)
					//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
					TDNR(transmission, phase, freq, dsc, num_iter)
				elseif (PopVal==2)
					//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
					TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
				endif
				
				//---Declare n and k and then set their scale---//
				Wave n,k
				SetScale/P x 0,(delta),"", n,k
				
				//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
				Duplicate/O $("n"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n:$("nr_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("n")
				
				Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n:$("kr_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:k:$("kr_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("k")
				
				Duplicate/O $("phase"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Phase:$("Thetar_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("phase")
					
			endfor
				
					//-----A bit of wave clean up-----//
					Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
		endfor
		
		for (i=0; i<NumTemps; i+=1)	
			SetDataFolder root:$Sample_Name + "_InField_Analysis"
					
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
			
			for (j=0; j<NumFields; j+=1)
				
				//---Name some strings for waves we'll need later---//
				fname1 = ("nr_" + Temps[i] + "_" + Fields[j] + "kG")
				fname2 = ("kr_" + Temps[i] + "_" + Fields[j] + "kG")
				fname3 = ("nr_" + RefTemp + "_" + Fields[j] + "kG")
				fname4 = ("kr_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
				fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname7 = ("Chi2r_" + Temps[i] + "_" + Fields[j] + "kG")
				fname8 = ("Chi1r_" +  Temps[i] + "_" + Fields[j] + "kG")
				
				fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
				
				SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$(RefTemp):RightHand:n
				Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n:$fname4
				
				SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:n
				
				//---Make some waves we'll need for the calculation---//
				Make/C/O/N=(numpnts($fname1)), $fname5
				Make/C/O/N=(numpnts($fname1)), $fname6
				Make/O/N=(numpnts($fname1)), $fname7
				Make/O/N=(numpnts($fname1)), $fname8
				Make/C/O/N=(numpnts($fname1)), $fname9
				
				//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
				Wave F1 = $fname1
				Wave F2 = $fname2
				Wave F3 = $fname3
				Wave F4 = $fname4
				Wave/C F5 = $fname5
				Wave/C F6 = $fname6
				Wave F7 = $fname7
				Wave F8 = $fname8
				Wave/C F9 = $fname9
			
				delta = DimDelta($fname1,0)
				
				SetScale/P x 0,(delta),"", $fname5
				SetScale/P x 0,(delta),"", $fname6
				SetScale/P x 0,(delta),"", $fname7
				SetScale/P x 0,(delta),"", $fname8
				SetScale/P x 0,(delta),"", $fname9
				
				F6 = F3 //+ sqrt(-1)*F4
				F5 = F1 + sqrt(-1)*F2
				
				F9 = (F5/F6)^2 - 1
				
				F8 = real(F9)
				F7 = imag(F9)
		
			endfor
				
				//---A loop which will do some wave clean up for us---//
				
				if (cmpstr(Temps[i], RefTemp) == 0)
					for (j=0; j<NumFields; j+=1)
					
						fname2 = ("kr_" + Temps[i] + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname7 = ("Chi2r_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1r_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
						
					endfor
					
				else
					for (j=0; j<NumFields; j+=1)
						fname2 = ("kr_" + Temps[i] + "_" + Fields[j] + "kG")
						fname3 = ("nr_" + RefTemp + "_" + Fields[j] + "kG")
						fname4 = ("kr_" + RefTemp + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
						fname7 = ("Chi2r_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1r_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:RightHand:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname3, $fname4, $fname5, $fname6, $fname7, $fname8, $fname9
					endfor
				endif
			endfor
		
			for (i=0; i<NumTemps; i+=1)
			
				//---Set the data folder to the main folder to find "Fields_Temps" waves
				SetDataFolder root:$Sample_Name + "_InField_Analysis"
			
				//---Creates folders for the complex conductivity, complex n, and phase---//
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:k
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Phase
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi1
				NewDataFolder/O root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi2
				
				fname = "Fields_" + Temps[i]
				Wave/T Fields = $fname
				NumFields = NumPnts($fname)
			
				//--- j runs over the number of Fields we have---//
				for (j=0; j<NumFields; j+=1)
			
					//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Transmissions
					index_calc("Tl_" + Temps[i] + "_" + Fields[j] + "kG")
					Wave transmission, phase, freq
					
					//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
					CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
					Wave W_coef
					phase -= W_coef[0]
					
					//-----Run the newton-raphson code to find the index-----//
					delta = DimDelta($("Tl_" + Temps[i] + "_" + Fields[j] + "kG"),0)
					
					//---Runs different version of the Newton-Raphson code depending on your selection to "Sample:" in the Analysis Palette"
					//----PopVal=1 corresponds to normal single crystal on an aperture, PopVal=2 corresponds to single crystal on a substate---//
					if (PopVal==1)
						//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
						TDNR(transmission, phase, freq, dsc, num_iter)
					elseif (PopVal==2)
						//---Advanced version of the code to include extra reflections and transmissions of a thick crystal on a substrate---//
						TDNR_Chris_nc(transmission, phase, freq, dsc, num_iter)
					endif
					
					//---Declare n and k and then set their scale---//
					Wave n,k
					SetScale/P x 0,(delta),"", n,k
					
					//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
					Duplicate/O $("n"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n:$("nl_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("n")
					
					Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n:$("kl_" + Temps[i] + "_" + Fields[j] + "kG")
					Duplicate/O $("k"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:k:$("kl_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("k")
					
					Duplicate/O $("phase"), root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Phase:$("Thetal_" + Temps[i] + "_" + Fields[j] + "kG")
					KillWaves $("phase")
						
				endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
			endfor
			
			for (i=0; i<NumTemps; i+=1)	
				SetDataFolder root:$Sample_Name + "_InField_Analysis"
						
				fname = "Fields_" + Temps[i]
				Wave/T Fields = $fname
				NumFields = NumPnts($fname)
				
				for (j=0; j<NumFields; j+=1)
					
					//---Name some strings for waves we'll need later---//
					fname1 = ("nl_" + Temps[i] + "_" + Fields[j] + "kG")
					fname2 = ("kl_" + Temps[i] + "_" + Fields[j] + "kG")
					fname3 = ("nl_" + RefTemp + "_" + Fields[j] + "kG")
					fname4 = ("kl_" + RefTemp + "_" + Fields[j] + "kG")
					
					fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
					fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
					
					fname7 = ("Chi2l_" + Temps[i] + "_" + Fields[j] + "kG")
					fname8 = ("Chi1l_" +  Temps[i] + "_" + Fields[j] + "kG")
					
					fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
					
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$(RefTemp):LeftHand:n
					Duplicate/O $fname3, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n:$fname3
					Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n:$fname4
					
					SetDataFolder root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:n
					
					//---Make some waves we'll need for the calculation---//
					Make/C/O/N=(numpnts($fname1)), $fname5
					Make/C/O/N=(numpnts($fname1)), $fname6
					Make/O/N=(numpnts($fname1)), $fname7
					Make/O/N=(numpnts($fname1)), $fname8
					Make/C/O/N=(numpnts($fname1)), $fname9
					
					//---Declare some dummy waves for the calculations of Complex Chi and Complex n---//
					Wave F1 = $fname1
					Wave F2 = $fname2
					Wave F3 = $fname3
					Wave F4 = $fname4
					Wave/C F5 = $fname5
					Wave/C F6 = $fname6
					Wave F7 = $fname7
					Wave F8 = $fname8
					Wave/C F9 = $fname9
				
					delta = DimDelta($fname1,0)
					
					SetScale/P x 0,(delta),"", $fname5
					SetScale/P x 0,(delta),"", $fname6
					SetScale/P x 0,(delta),"", $fname7
					SetScale/P x 0,(delta),"", $fname8
					SetScale/P x 0,(delta),"", $fname9
					
					F6 = F3 //+ sqrt(-1)*F4
					F5 = F1 + sqrt(-1)*F2
					
					F9 = (F5/F6)^2 - 1
					
					F8 = real(F9)
					F7 = imag(F9)
			
				endfor
					
					//---A loop which will do some wave clean up for us---//
					
					if (cmpstr(Temps[i], RefTemp) == 0)
						for (j=0; j<NumFields; j+=1)
						
							fname2 = ("kl_" + Temps[i] + "_" + Fields[j] + "kG")
							fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
							fname7 = ("Chi2l_" + Temps[i] + "_" + Fields[j] + "kG")
							fname8 = ("Chi1l_" +  Temps[i] + "_" + Fields[j] + "kG")
							fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
							
							Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi2:$fname7
							Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi1:$fname8
							
							Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
							
						endfor
						
					else
						for (j=0; j<NumFields; j+=1)
							fname2 = ("kl_" + Temps[i] + "_" + Fields[j] + "kG")
							fname3 = ("nl_" + RefTemp + "_" + Fields[j] + "kG")
							fname4 = ("kl_" + RefTemp + "_" + Fields[j] + "kG")
							fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
							fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
							fname7 = ("Chi2l_" + Temps[i] + "_" + Fields[j] + "kG")
							fname8 = ("Chi1l_" +  Temps[i] + "_" + Fields[j] + "kG")
							fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
							
							Duplicate/O $fname7, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi2:$fname7
							Duplicate/O $fname8, root:$Sample_Name + "_InField_Analysis":Circular_Basis:$Temps[i]:LeftHand:Chi1:$fname8
							
							Killwaves/Z $fname2, $fname3, $fname4, $fname5, $fname6, $fname7, $fname8, $fname9
						endfor
					endif
			endfor
		endif
	endif
endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the Faraday Rotation for either the BLC or the Magnet.  Only works when the rotator option is chosen.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function FaradayRotation()

//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----Dummy variables-----//
	String fname, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9, fname10;
	variable i,k,j;
	
	ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
		if (V_Value ==1)
		
			//---This is the case that the "BLC" option is selected on the front panel---//
			ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
			if (V_Value == 1)
						
				//---Sets the data folder to the Analysis folder---//
				SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
				
				//---recognizes the waves that we'll need while running this function---//
				Wave/T Temps
				Wave num_samp_scans
				Variable NumTemps = numpnts(Temps)
				
				//---Creates folders for the transmissions and the magntiude of the transmissoins---//
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation'
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':'Faraday Angle'
				NewDataFolder/O/S root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':Ellipticity 
				
				for(i=0; i<NumTemps; i+=1)
					
					//---These strings hold the linear basis transmissions---///
					fname = ("Txx_" + Sample_Name + "_" + Temps[i])
					fname2 = ("Tyx_" + Sample_Name + "_" + Temps[i])
					
					//---These strings hold the Faraday Rotation Wave---//
					fname3 = ("FR_" + Sample_Name + "_" + Temps[i])
					fname4 = ("FR1_" + Sample_Name + "_" + Temps[i])
					fname5 = ("FR2_" + Sample_Name + "_" + Temps[i])

					//---First go to the in phase transmissions and Duplicate Txx---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:Transmissions
					Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':$fname
					
					//---Then go to the our of phase transmissions and Duplicate Tyx---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":Out_Of_Phase:Transmissions
					Duplicate/O $fname2, root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':$fname2
					
					//---Duplicate some wave and rename them as the ciruclar basis transmissions---//
					SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation'
					Duplicate/O $fname2, $fname3
					Duplicate/O $fname2, $fname4
					Duplicate/O $fname2, $fname5
					
					Redimension/R $fname4
					Redimension/R $fname5
					
					//---Declare some dummy waves for the calculation---//
					Wave/C DummyWave = $fname
					Wave/C DummyWave2 = $fname2
					Wave/C DummyWave3 = $fname3
					Wave DummyWave4 = $fname4
					Wave DummyWave5 = $fname5
					
					//---Calculate the Faraday Rotation---//
					DummyWave3 = atan(DummyWave2 / DummyWave)
					
					DummyWave4 = Real(DummyWave3)
					DummyWave5 = Imag(DummyWave3)
					
					Duplicate/O $fname4,  root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':'Faraday Angle':$fname4
					Duplicate/O $fname5,  root:$Sample_Name + "_ZeroField_Analysis":'Faraday Rotation':Ellipticity:$fname5
					
					KillWaves/Z $fname, $fname2, $fname3, $fname4, $fname5
					
				endfor	
			endif
				
			//---This is the case that we're using the magnet---//
			ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
			if (V_Value == 1)
			
				//---Define some variables---//
				variable NumFields, p;
			
				SetDataFolder root:$Sample_Name + "_InField_Analysis":
				Wave/T Temps
				NumTemps = NumPnts(Temps)
				
				for (i=0; i<NumTemps; i+=1)
				
					//---Set the data folder to the main folder to find "Fields_Temps" waves
					SetDataFolder root:$Sample_Name + "_InField_Analysis":
						
					//---Creates folders for the transmissions and the magntiude of the transmissoins---//
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Faraday Rotation'
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:'Faraday Angle'
					NewDataFolder/O root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:Ellipticity
					
					fname = "Fields_" + Temps[i]
					Wave/T Fields = $fname
					NumFields = NumPnts($fname)
				
					//--- j runs over the number of Fields we have---//
					for (j=0; j<NumFields; j+=1)
					
						//---These strings hold the linear basis transmissions---//
						fname = ("Txx_" + Temps[i] + "_" + Fields[j] + "kG")
						fname2 = ("Tyx_" + Temps[i] + "_" + Fields[j] + "kG")
						
						//---These strings hold the circular basis transmissions---//
						fname3 = ("FR_" + Temps[i] + "_" + Fields[j] + "kG")
						fname4 = ("FR1_" + Temps[i] + "_" + Fields[j] + "kG")
						fname5 = ("FR2_" + Temps[i] + "_" + Fields[j] + "kG")
	
						//---First go to the in phase transmissions and Duplicate/O Txx---//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":In_Phase:$Temps[i]:Transmissions
						Duplicate/O $fname, root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:$fname
						
						//---Then go to the our of phase transmissions and Duplicate/O Tyx---//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":Out_Of_Phase:$Temps[i]:Transmissions
						Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:$fname2
						
						//---Duplicate/O some wave and rename them as the ciruclar basis transmissions---//
						SetDataFolder root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]
						Duplicate/O $fname2, $fname3
						Duplicate/O $fname2, $fname4
						Duplicate/O $fname2, $fname5
						
						Redimension/R $fname4, $fname5
						
						//---Declare some dummy waves for the calculation---//
						Wave/C DummyWave = $fname
						Wave/C DummyWave2 = $fname2
						Wave/C DummyWave3 = $fname3
						Wave DummyWave4 = $fname4
						Wave DummyWave5 = $fname5
						
						//---Calculate Faraday Rotation---//
						DummyWave3 = Atan(DummyWave2 / DummyWave)
						
						DummyWave4 = Real(DummyWave3)
						DummyWave5 = Imag(DummyWave3)
						
						Duplicate/O $fname4, root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:'Faraday Angle':$fname4	
						Duplicate/O $fname5, root:$Sample_Name + "_InField_Analysis":'Faraday Rotation':$Temps[i]:Ellipticity:$fname5
						
						//---Delete the linear transmissions---//
						KillWaves/Z $fname, $fname2, $fname3, $fname4, $fname5
			
					endfor	
				endfor		
			endif
	else
		Print "You must use the rotator to calculate Faraday rotation!"
	endif
End


//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function is used to create folders with data taken at a constant field.  Folder1 can be either "In_Phase" or "Out_Of_Phase", Folder2 is what you want to call the folder in the constant field folder.
// For instance "Transmission" or "Sigma1" etc.  Folder 3 is the folder you want to import the data from, this will almost always be the same as Folder 2 except for the case of Transmissions.  WaveTitle is the name of the
// wave which again will often be the same as Folders 2 and 3.  However it can be different.  For instance, for the phase Folders 2 and 3 will be "Phase" while WaveTitle will be "Theta"
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function ConstantField(Folder1, Folder2, Folder3, Wavetitle)
//---Declare our parameters---//
	String Folder1, Folder2, Folder3, Wavetitle
	
//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
	
//-----number of Temps that were measured at------//
	NVAR num_temps = root:Packages:NickLTHzExtras:num_temps

//-----Dummy variables-----//
	String fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9, fname10, OtherFolder1, OtherFolder2, WaveEnding, OtherWaveEnding;
	variable i,k,j;
		
//---Only generate constant field folders if we're running the magnet---//
ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
if (V_Value ==1)
	
	//---This if statement sets the values of "OtherFolder" 1 and 2 according you what you put into Folder 1 and 2---//
	//---Essentially makes sure the in phase and out of phase folders are creates, to be honest, it's pretty redundant---//
	if (cmpstr(Folder1, "In_Phase") == 0)
	
		OtherFolder1 = "Out_of_Phase"
		OtherFolder2 = Folder2
		WaveEnding = "xx"
		OtherWaveEnding = "yx"
		
	elseif (cmpstr(Folder1, "Out_of_Phase") == 0)
	
		OtherFolder1 = "In_Phase"
		OtherFolder2 = Folder2
		WaveEnding = "yx"
		OtherWaveEnding = "xx"
		
	elseif (cmpstr(Folder1, "Circular_Basis") == 0)
	
		if (cmpstr(Folder2, "RightHand") == 0)
			OtherFolder1 = Folder1
			OtherFolder2 = "LeftHand"
			WaveEnding = "r"
			OtherWaveEnding = "l"
			
		elseif (cmpstr(Folder2, "LeftHand") == 0)
			OtherFolder1 = Folder1
			OtherFolder2 = "RightHand"
			WaveEnding = "l"
			OtherWaveEnding = "r"
		endif
	endif
	
	//---Define some variables---//
	variable NumFields, p;
	NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(Folder1):Constant_Field
	NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(Folder1):Constant_Field:$(Folder2)
	
	//---Checks to see if we're using the rotator---//
	ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
	if (V_Value ==1)
		
		NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(OtherFolder1):Constant_Field
		NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(OtherFolder1):Constant_Field:$(OtherFolder2)
		
	endif

	SetDataFolder root:$Sample_Name + "_InField_Analysis":
	Wave/T Temps
	Variable NumTemps = NumPnts(Temps)
	
	for (i=0; i<NumTemps; i+=1)
	
		SetDataFolder root:$Sample_Name + "_InField_Analysis":
		String fname = "Fields_" + Temps[i]
		Wave/T Fields = $fname
		NumFields = NumPnts($fname)
	
		//--- j runs over the number of Fields we have---//
		for (j=0; j<NumFields; j+=1)

					//---Create a data folder for each field that's there---//
					NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(Folder1):Constant_Field:$(Folder2):$Fields[j]+"_kG"
					
					if (cmpstr(Folder1, "Faraday Rotation") == 0)
						fname1 = WaveTitle  + "_" + Temps[i] + "_" + Fields[j] + "kG"
					else
						fname1 = WaveTitle + WaveEnding + "_" + Temps[i] + "_" + Fields[j] + "kG"
					endif
		
					// Grabs the In Phase $(Folder3)
					SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$Temps[i]:$(Folder3)
					Duplicate/O $fname1, root:$Sample_Name + "_InField_Analysis":$(Folder1):Constant_Field:$(Folder2):$Fields[j]+"_kG":$fname1
					
					//---Checks to see if we're using the rotator---//
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value ==1)
		
						NewDataFolder/O/S root:$Sample_Name + "_InField_Analysis":$(OtherFolder1):Constant_Field:$(OtherFolder2):$Fields[j]+"_kG"
						
						if (cmpstr(Folder1, "Faraday Rotation") == 0)
							fname2 = WaveTitle + "_" + Temps[i] + "_" + Fields[j] + "kG"
						else	
							fname2 = WaveTitle + OtherWaveEnding + "_" + Temps[i] + "_" + Fields[j] + "kG"
						endif
						
						// Grabs the Out of Phase $(Folder3)
						SetDataFolder root:$Sample_Name + "_InField_Analysis":$(OtherFolder1):$Temps[i]:$(Folder3)
						Duplicate/O $fname2, root:$Sample_Name + "_InField_Analysis":$(OtherFolder1):Constant_Field:$(OtherFolder2):$Fields[j]+"_kG":$fname2
						
					endif
		
		endfor
	endfor
endif
//endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function uses the previous function that makes the constant field folders but checks to see which type of analysis has been done so far and generates constant field plots accordingly.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function GenerateConstantField()
	
//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder

//-----Dummy variables-----//
	String fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9, fname10, OtherFolder1, OtherFolder2, WaveEnding, OtherWaveEnding;
	variable i,k,j;
	
	SetDataFolder root:$Sample_Name + "_InField_Analysis":
	Wave/T Temps
	String FirstTemp = Temps[0]
	
		//---Always generate the constant field folder for the Transmissions---//
		ConstantField("In_Phase", "Transmissions", "Mag_Transmissions", "Mag_T")
		
		if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:In_Phase:'" + FirstTemp + "':n:") == 1)
			ConstantField("In_Phase", "n", "n", "n")
			ConstantField("In_Phase", "k", "k", "k")
			ConstantField("In_Phase", "Phase", "Phase", "Theta")
		endif
		
			if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:In_Phase:'" + FirstTemp + "':Sigma1:") == 1)
				ConstantField("In_Phase", "Sigma1", "Sigma1", "Sigma1")
				ConstantField("In_Phase", "Sigma2", "Sigma2", "Sigma2")
			endif
		
			if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:In_Phase:'" + FirstTemp + "':Conductance1:") == 1)
				ConstantField("In_Phase", "Conductance1", "Conductance1", "G1")
				ConstantField("In_Phase", "Conductance2", "Conductance2", "G2")
			endif
		
			if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:In_Phase:'" + FirstTemp + "':Chi1:") == 1)
				ConstantField("In_Phase", "Chi1", "Chi1", "Chi1")
				ConstantField("In_Phase", "Chi2", "Chi2", "Chi2")
			endif

			if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:'Faraday Rotation'") == 1)
				ConstantField("Faraday Rotation", "Faraday Angle", "Faraday Angle", "FR1")
				ConstantField("Faraday Rotation", "Ellipticity", "Ellipticity", "FR2")
			endif
			
//			if (DataFolderExists("root:" + Sample_Name + "_InField_Analysis:Circular_Basis") == 1)
//				ConstantField("Circular_Basis", "RightHand", "Mag_Transmissions", "Mag_T")
//				ConstantField("Faraday Rotation", "Ellipticity", "Ellipticity", "FR2")
//			endif

End


Function FreqCutsSCZeroField()
//---This is a function that takes the frequency cuts if you've selected one of the thick crystal opetions to "Sample:" in the Analysis Palette"

	//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
	
	//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
	//-----Starting Frequency for the cuts------//
	NVAR Freq_Low_Limit = root:Packages:NickLTHzExtras:Freq_Low_Limit

	//-----Ending Frequency for the cuts-----//
	NVAR Freq_High_Limit = root:Packages:NickLTHzExtras:Freq_High_Limit

	//-----Step size for the frequency cuts------//
	NVAR Cut_Resolution = root:Packages:NickLTHzExtras:Cut_Resolution

	//---Sets the data folder to the Analysis folder just in case---//
	SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
	
	//---recognizes the wave that we'll need while running this function---//
	Wave/T Temps
	Wave num_samp_scans
	Variable NumTemps = numpnts(Temps)
	
	//-----Dummy variables-----//
	String fname, fname2;
	
	//---Declare some variables---//
	Variable i, j, k, Delta;
	Variable cut_iter, current_freq, num_cut_wave_types
	
	if (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:In_Phase:Sigma1:") == 1)
		Print "oh"

		//----These perform constant frequency cuts as a function of temperature-----//
		//----You can grab the magnitude, real, or imaginary parts of any wave series by adding "mag", "real", or "imag" to CutWaveType
		num_cut_wave_types = 5
		Make/O/T/N=(num_cut_wave_types) CutWaveBaseNames
		Make/O/T/N=(num_cut_wave_types) FolderNames
		Make/O/T/N=(num_cut_wave_types) CutWaveType //Can be magnitude ("mag"), real part ("real"), imaginary part ("imag")
		//-----If the wave is only real, mag or real will work, but not imag-----//
		
		CutWaveBaseNames[0] = "T";						CutWaveType[0] = "mag";				FolderNames[0] = "Transmissions"
		CutWaveBaseNames[1] = "sigma1"; 				CutWaveType[1] = "mag";				FolderNames[1] = "Sigma1"
		CutWaveBaseNames[2] = "sigma2"; 				CutWaveType[2] = "mag";				FolderNames[2] = "Sigma2"
		CutWaveBaseNames[3] = "n"; 						CutWaveType[3] = "mag";				FolderNames[3] = "n"
		CutWaveBaseNames[4] = "k";						CutWaveType[4] = "mag";				FolderNames[4] = "k"
	
	elseif (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:In_Phase:Conductance1:") == 1)
		Print "my"
	
		//----These perform constant frequency cuts as a function of temperature-----//
		//----You can grab the magnitude, real, or imaginary parts of any wave series by adding "mag", "real", or "imag" to CutWaveType
		num_cut_wave_types = 3
		Make/O/T/N=(num_cut_wave_types) CutWaveBaseNames
		Make/O/T/N=(num_cut_wave_types) FolderNames
		Make/O/T/N=(num_cut_wave_types) CutWaveType //Can be magnitude ("mag"), real part ("real"), imaginary part ("imag")
		//-----If the wave is only real, mag or real will work, but not imag-----//
		
		CutWaveBaseNames[0] = "T";						CutWaveType[0] = "mag";				FolderNames[0] = "Transmissions"
		CutWaveBaseNames[1] = "G1"; 					CutWaveType[1] = "mag";				FolderNames[1] = "Conductance1"
		CutWaveBaseNames[2] = "G2"; 					CutWaveType[2] = "mag";				FolderNames[2] = "Conductance2"
		
	elseif (DataFolderExists("root:" + Sample_Name + "_ZeroField_Analysis:In_Phase:Chi1:") == 1)
		Print "god"
	
		//----These perform constant frequency cuts as a function of temperature-----//
		//----You can grab the magnitude, real, or imaginary parts of any wave series by adding "mag", "real", or "imag" to CutWaveType
		num_cut_wave_types = 5
		Make/O/T/N=(num_cut_wave_types) CutWaveBaseNames
		Make/O/T/N=(num_cut_wave_types) FolderNames
		Make/O/T/N=(num_cut_wave_types) CutWaveType //Can be magnitude ("mag"), real part ("real"), imaginary part ("imag")
		//-----If the wave is only real, mag or real will work, but not imag-----//
		
		CutWaveBaseNames[0] = "T";						CutWaveType[0] = "mag";				FolderNames[0] = "Transmissions"
		CutWaveBaseNames[1] = "Chi1"; 					CutWaveType[1] = "mag";				FolderNames[1] = "Chi1"
		CutWaveBaseNames[2] = "Chi2"; 					CutWaveType[2] = "mag";				FolderNames[2] = "Chi2"
		CutWaveBaseNames[3] = "n"; 						CutWaveType[3] = "mag";				FolderNames[3] = "n"
		CutWaveBaseNames[4] = "k";						CutWaveType[4] = "mag";				FolderNames[4] = "k"
	
	endif
	
	//---Duplicate Sigma1, Sigma2, n, k and move them to the main folder.
	for(cut_iter = 0; cut_iter < num_cut_wave_types; cut_iter += 1)
		for(j=0; j<numpnts(Temps); j+=1)
		
			fname = CutWaveBaseNames[cut_iter] + "xx_" + Sample_Name + "_" + Temps[j]
			fname2 = FolderNames[cut_iter]
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase:$fname2 
			Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:$fname
			
		endfor
	endfor
	
	for(cut_iter = 0; cut_iter < num_cut_wave_types; cut_iter += 1)
	
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":In_Phase
	
		//-----Create the data folder to hold the cuts----//
		String data_folder_name = "Freq_Cuts_" + CutWaveBaseNames[cut_iter] //+ "_" + CutWaveType[cut_iter]
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":In_Phase:$data_folder_name

		String base_file_name = CutWaveBaseNames[cut_iter] + "xx_" + Sample_Name
		String choose_type = CutWaveType[cut_iter]
		Variable num_cuts = (freq_high_limit - freq_low_limit)/cut_resolution + 1
		
		//-----This section duplicates the desired part of the current wave, but as a real wave----//
		//-----We do this because we want to be able to interpolate to get the desired frequency point----//
		//-----and interp() only works on non-complex waves-----//
		Make/O/N = (numpnts(Temps)) Temps_Num;
		for(i = 0; i < numpnts(Temps); i += 1)
			fname = base_file_name + "_" + Temps[i]
			Wave current_trans = $fname
			Make/O/N=(numpnts(current_trans)), $(base_file_name + "_Dummy_" + Temps[i])
			Wave tdummy =  $(base_file_name + "_Dummy_" + Temps[i]) 
			if(cmpstr(choose_type,"mag") == 0)
				tdummy = cabs(current_trans);
			elseif(cmpstr(choose_type,"real") == 0)
				tdummy = real(current_trans);
			elseif(cmpstr(choose_type,"imag") == 0)
				tdummy = imag(current_trans);
			endif
			delta = DimDelta(current_trans,0); SetScale/P x 0,(delta),"", tdummy
			Temps_Num[i] = str2num(Temps[i])
		endfor
		
		//-----interp() also needs the x wave, it wont work just knowing the delta of the given y wave, so I create this-----//
		fname = base_file_name + "_" + Temps[0]
		delta = DimDelta($fname,0)
		Make/O/N=(numpnts($fname)) frequency
		for(i = 0; i < numpnts(frequency); i += 1)
			frequency[i] = delta * i
		endfor
		
		//-----Create the wave and store the temperature dependent cut values at each frequency-----//
		for(i = 0; i < num_cuts; i += 1)
			current_freq = cut_resolution * i + freq_low_limit   //Current cut frequency
			//-----Creates the wave to hold the temperature dep. data at the current frequency-----//
			Make/O/N=(numpnts(Temps)) $(base_file_name + "_" + num2str(current_freq) + "THz")
			Wave dummy_fcut = $(base_file_name + "_" + num2str(current_freq) + "THz")
			//-----At each frequency run through all the Temps and store the value cut we are interested in such as transmission, sigma, etc.-----//
			for(j = 0; j<numpnts(Temps); j += 1)
				Wave current_base_wave = $(base_file_name + "_Dummy_" + Temps[j])
				dummy_fcut[j] = interp(current_freq, frequency, current_base_wave)  //We have to interpolate since we force evenly spaced values of x we decide on
			endfor
		endfor
		
		//-----Move all these generated cuts at each frequency into the appropriate folder-----//
		Duplicate/O Temps_Num, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:$(data_folder_name):Temps
		for(i = 0; i < num_cuts; i += 1)
			current_freq = cut_resolution * i + freq_low_limit
			fname = base_file_name + "_" + num2str(current_freq) + "THz"
			Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":In_Phase:$(data_folder_name):$(fname) //+ "_" + CutWaveType[cut_iter])
			KillWaves $fname
		endfor
		
		//-----Kill all the dummy real waves we had to create to get interp() to work above-----//
		for(i = 0; i < numpnts(Temps); i += 1)
			KillWaves  $(base_file_name + "_Dummy_" + Temps[i])
		endfor
	endfor
	
	//---Now it's just a bit of wave clean up.  Let's kill all the waves that we put into the main folder, we don't need them there anymore.
	SetDataFolder  root:$Sample_Name + "_ZeroField_Analysis":In_Phase
	for(i=0; i<num_cut_wave_types; i+=1)
		for(j=0; j<numpnts(Temps); j+=1)
				fname = CutWaveBaseNames[i] + "xx_" + Sample_Name + "_" +  Temps[j]
				KillWaves/Z $fname
		endfor
	endfor
	
	//---Set the data folder back the main folder for the next function---//
	SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
	
	//---Lets kill the wave associated with taking the cuts.
	KillWaves/Z CutWaveBaseNames, CutWaveType, Temps_Num, FolderNames

	
End

Function FreqCutsTFZeroField()
//---This function performs the frequency cuts in the case of a thin film sample

		//-----Folder where data is loaded from-----//
		SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
		
		//-----String where you enter in the name of the sample and aperature-----//
		SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
		SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
		
		//-----number of Temps that were measured at------//
		NVAR num_temps = root:Packages:NickLTHzExtras:num_temps
		
		//-----Starting Frequency for the cuts------//
		NVAR Freq_Low_Limit = root:Packages:NickLTHzExtras:Freq_Low_Limit
	
		//-----Ending Frequency for the cuts-----//
		NVAR Freq_High_Limit = root:Packages:NickLTHzExtras:Freq_High_Limit
	
		//-----Step size for the frequency cuts------//
		NVAR Cut_Resolution = root:Packages:NickLTHzExtras:Cut_Resolution

	//---Sets the data folder to the Analysis folder just in case---//
	SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
	
	//---recognizes the wave that we'll need while running this function---//
	Wave/T Temps
	Wave num_samp_scans
	num_temps = numpnts(Temps)
	
	//-----Dummy variables-----//
	String fname, fname2;
	
	Variable i, j, k, Delta;

	//----These perform constant frequency cuts as a function of temperature-----//
	//----You can grab the magnitude, real, or imaginary parts of any wave series by adding "mag", "real", or "imag" to CutWaveType
	Variable num_cut_wave_Types = 3
	Make/O/T/N=(num_cut_wave_Types) CutWaveBaseNames
	Make/O/T/N=(num_cut_wave_Types) CutWaveType //Can be magnitude ("mag"), real part ("real"), imaginary part ("imag")
	//-----If the wave is only real, mag or real will work, but not imag-----//
	
	CutWaveBaseNames[0] = "T_"+ Sample_Name;	 			CutWaveType[0] = "mag"
	CutWaveBaseNames[1] = "Conductance1"; 					CutWaveType[1] = "mag"
	CutWaveBaseNames[2] = "Conductance2"; 					CutWaveType[2] = "mag"

	
	Variable cut_iter, current_freq
	
	//---In order for the code to work correctly I'm going to copy all the waves that we want to take cuts of from their folders and duplicate them 
	//---in the main folder.  The transmissions are a little tricky because the name of the folder is "transmissions" but the files are "T_Sample_Name..".
	//---So the first loop duplicates the transmissions while the second loop duplicated everything else.
	
	//---Duplicate the transmissions and move them to the main folder
	for(j=0; j<num_temps; j+=1)
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
		fname = CutWaveBaseNames[0] + "_" + Temps[j]
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":'Transmissions'
		Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$fname
		
	endfor
	
	//---Duplicate Sigma1, Sigma2, n, k and move them to the main folder.
	for(cut_iter = 1; cut_iter < num_cut_wave_Types; cut_iter += 1)
		for(j=0; j<numpnts(Temps); j+=1)
		
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
			fname = CutWaveBaseNames[cut_iter] + "_" + Temps[j]
			fname2 = CutWaveBaseNames[cut_iter]
			SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$fname2 
			Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$fname
			
		endfor
	endfor
	
	for(cut_iter = 0; cut_iter < num_cut_wave_Types; cut_iter += 1)
	
		SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
	
		//-----Create the data folder to hold the cuts----//
		String data_folder_name = "Freq_Cuts_" + CutWaveBaseNames[cut_iter]// + "_" + CutWaveType[cut_iter]
		NewDataFolder/O root:$Sample_Name + "_ZeroField_Analysis":$data_folder_name
		
		String base_file_name = CutWaveBaseNames[cut_iter] + "_"
		String choose_type = CutWaveType[cut_iter]
		Variable num_cuts = (freq_high_limit - freq_low_limit)/cut_resolution + 1
		
		//-----This section duplicates the desired part of the current wave, but as a real wave----//
		//-----We do this because we want to be able to interpolate to get the desired frequency point----//
		//-----and interp() only works on non-complex waves-----//
		Make/O/N = (num_temps) Temps_Num;
		for(i = 0; i < num_temps; i += 1)
			fname = base_file_name + Temps[i]
			Wave current_trans = $fname
			Make/O/N=(numpnts(current_trans)), $(base_file_name + "Dummy_" + Temps[i])
			Wave tdummy =  $(base_file_name + "Dummy_" + Temps[i])
			if(cmpstr(choose_type,"mag") == 0)
				tdummy = cabs(current_trans);
			elseif(cmpstr(choose_type,"real") == 0)
				tdummy = real(current_trans);
			elseif(cmpstr(choose_type,"imag") == 0)
				tdummy = imag(current_trans);
			endif
			delta = DimDelta(current_trans,0); SetScale/P x 0,(delta),"", tdummy
			Temps_Num[i] = str2num(Temps[i])
		endfor
		
		//-----interp() also needs the x wave, it wont work just knowing the delta of the given y wave, so I create this-----//
		fname = base_file_name + Temps[0]
		delta = DimDelta($fname,0)
		Make/O/N=(numpnts($fname)) frequency
		for(i = 0; i < numpnts(frequency); i += 1)
			frequency[i] = delta * i
		endfor
		
		//-----Create the wave and store the temperature dependent cut values at each frequency-----//
		for(i = 0; i < num_cuts; i += 1)
			current_freq = cut_resolution * i + freq_low_limit   //Current cut frequency
			//-----Creates the wave to hold the temperature dep. data at the current frequency-----//
			Make/O/N=(num_temps) $(base_file_name + num2str(current_freq) + "THz")
			Wave dummy_fcut = $(base_file_name + num2str(current_freq) + "THz")
			//-----At each frequency run through all the Temps and store the value cut we are interested in such as transmission, sigma, etc.-----//
			for(j = 0; j< num_temps; j += 1)
				Wave current_base_wave = $(base_file_name + "Dummy_" + Temps[j])
				dummy_fcut[j] = interp(current_freq, frequency, current_base_wave)  //We have to interpolate since we force evenly spaced values of x we decide on
			endfor
		endfor
		
		//-----Move all these generated cuts at each frequency into the appropriate folder-----//
		Duplicate/O Temps_Num, root:$Sample_Name + "_ZeroField_Analysis":$(data_folder_name):Temps
		for(i = 0; i < num_cuts; i += 1)
			current_freq = cut_resolution * i + freq_low_limit
			fname = base_file_name + num2str(current_freq) + "THz"
			Duplicate/O $fname, root:$Sample_Name + "_ZeroField_Analysis":$(data_folder_name):$(fname)// + "_" + CutWaveType[cut_iter])
			KillWaves $fname
		endfor
		
		//-----Kill all the dummy real waves we had to create to get interp() to work above-----//
		for(i = 0; i < num_temps; i += 1)
			KillWaves  $(base_file_name + "Dummy_" + Temps[i])
		endfor
	endfor
	
	//---Now it's just a bit of wave clean up.  Let's kill all the waves that we put into the main folder, we don't need them there anymore.
	for(i=0; i<num_cut_wave_Types; i+=1)
		for(j=0; j<numpnts(Temps); j+=1)
				fname = CutWaveBaseNames[i] + "_" +  Temps[j]
				KillWaves/Z $fname
		endfor
	endfor
	
	//---Lets kill the wave associated with taking the cuts.
	KillWaves/Z CutWaveBaseNames, CutWaveType
	
	//---Set the data folder back the main folder for the next function---//
	SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
	
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function takes in as parameters the folder and the general wave name and then populates a graph of the data.  Should work for anything that's calculated via this panel in both the BLC or Magnet with or without 
// the rotator.  Q is a parameter which permutates the possible paths for the data folder.  Folder1 refers to the "In_Phase" or "Out_Of_Phase" part of the wave path.  Folder2 is the next folder down.  
//You can trick the code into running without the rotator by putting the second folder in for folder1 and  then having nothing in folder2.  For instance MakePlot("Transmissions", "", "T") would run without the rotator.
// I added an additional folder, folder3, so as to make the code more versitile.  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function MakePlot(Q, Folder1, Folder2, Folder3, WaveTitle)
	//---Defines the number of folders we want to have in the title---//
	Variable Q;
	
	//---Define the folder and the wave that we want on each graph---//
	String WaveTItle, Folder1, Folder2, Folder3;
	
	//-----Folder where data is loaded from-----//
	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
	
	//-----String where you enter in the name of the sample and aperature-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
	
	//---Define some dummy variables---//
	String fname, fname2, fname3, fname4;
	variable i,j,k;
	
	//---Only create the graphs if the check box on the front panel is selected---//
		ControlInfo/W=NickL_THzAnalysis_Palette DoGraphs_CheckBox
			if (V_Value == 1)
	
			//---This is the case that the "BLC" option is selected on the front panel---//
				ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
					if (V_Value == 1)
	
						//---recognizes the waves that we'll need while running this function---//
						SetDataFolder root:$Sample_Name + "_ZeroField_Analysis"
						Wave/T Temps
						Wave num_samp_scans
						Variable NumTemps = numpnts(Temps)
						
						if (Q==1)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$(Folder1)
						elseif (Q==2)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):$(Folder2)
						elseif (Q==3)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_ZeroField_Analysis":$(Folder1):$(Folder2):$(Folder3)
						else
							abort
						endif
							
							//---Lets kill the window so we can recreate it, doesn't throw an error if the window doesn't exist yet---//
							DoWindow/K $(WaveTitle + "_vs_Frequency")
							//---Now let's create the window so we have make the graph---//
							Display/N= $(WaveTitle + "_vs_Frequency")
							
							//---Append the traces to the graph--//
							for(i=0; i<numpnts(Temps); i+=1)
								fname = (WaveTitle + "_" +  Sample_Name + "_" + Temps[i])
								AppendToGraph/W =  $(WaveTitle + "_vs_Frequency") $fname
							endfor
							
							//---Make it look like we want it to---///
							ModifyGraph/Z cmplxMode=3, lsize=2, mirror=2, fStyle=1, fSize=14, axThick=2, nticks=10;
							Setaxis/Z/W= $(WaveTitle + "_vs_Frequency") bottom 0.2, 2.25;
							Setaxis/Z/W= $(WaveTitle + "_vs_Frequency")/A=2 left; 
							Legend/C/N=text0/F=0/H={30,1,10}/A=MC
							Label bottom "\\Z16\\f01Frequency (THz)"
							Label left "\\Z16\\f01" + Folder1 + " " + Folder2
								
							KBColorizeTraces#KBColorTablePopMenuProc(WaveTitle +"_vs_Frequency",0,"BlueRedGreen")	
					endif
		
		//---This is the case that the "Magnet" option is selected on the front panel---//
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
			if (V_Value == 1)
			
				variable p=0, Num_TempFields;
				
				//--- This loop tells us how many temperatures were measured at in the Magnet system---//
				for (k=1; k<=100; k+=1)
					SetDataFolder root:$Sample_Name + "_InField_Analysis":
					fname = "TempFields_" + num2str(k)
					Wave/T TFNow = $fname
					p += WaveExists(TFNow)
				endfor
					Num_TempFields = p
					
					for (k=1; k<=Num_TempFields; k+=1)

						SetDataFolder root:$Sample_Name + "_InField_Analysis":
					
						String CurrentTF = "TempFields_" + num2str(k)
						Wave/T TFNow = $CurrentTF
						
						Variable NumFields = numpnts(TFNow) 
						
						if (Q==1)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":$TFNow[0]:$(Folder1)
						elseif (Q==2)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$TFNow[0]:$(Folder2)
						elseif (Q==3)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$TFNow[0]:$(Folder2):$(Folder3)
						elseif (Q==4)
							//---Sets the current folder---//
							SetDataFolder root:$Sample_Name + "_InField_Analysis":$(Folder1):$(Folder2):$TFNow[0]:$(Folder3)
						else
							abort
						endif
						
						//---Lets kill the window so we can recreate it, doesn't throw an error if the window doesn't exist yet---//
						DoWindow/K $(WaveTitle + "_vs_Frequency_" + TFNow[0])
						//---Now let's create the window so we have make the graph---//
						Display/N=$(WaveTitle + "_vs_Frequency_" + TFNow[0])

						//---Loop that runs over fields---//
						for( j = 1; j < NumFields; j += 1)

							//---Append the traces to the graph--//
								fname = (WaveTitle + "_" + TFNow[0] + "_" + TFNow[j] + "kG")
								AppendToGraph/W = $(WaveTitle + "_vs_Frequency_" + TFNow[0]) $fname					
							
								//---Make it look like we want it to---///
								ModifyGraph/Z cmplxMode=3, lsize=2, mirror=2, fStyle=1, fSize=14, axThick=2, nticks=10;
								Setaxis/Z/W= $(WaveTitle + "_vs_Frequency_" + TFNow[0]) bottom 0.2, 2.25;
								Setaxis/Z/W= $(WaveTitle + "_vs_Frequency_" + TFNow[0])/A=2 left; 
								Legend/C/N=text0/F=0/H={30,1,10}/A=MC
								Label bottom "\\Z16\\f01Frequency (THz)"
								Label left "\\Z16\\f01" + Folder1 + " " + Folder2
						
								//---Colorize the trances so they look very nice---//
								KBColorizeTraces#KBColorTablePopMenuProc(WaveTitle + "_vs_Frequency_" + TFNow[0],0,"BlueRedGreen")
						endfor
					endfor	
			endif
	endif
End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function makes a plot of all the data that's contained in a given folder.  Simply select that folder as the default folder in the data browser and it will upload all the data to a plot.  Can add offset as well. 
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function GenerateGraph()

	//-----OSA is the offset amount-----//
	NVAR OSA = root:Packages:NickLTHzExtras:OSA;
	
	//---Define some dummy variables---//
	String fname, fname2, fname3, fname4, WaveTitle, AxisTitle, WaveEnding, Phase
	variable i,j,k;
	
	String FullDataFolder
	FullDataFolder = GetDataFolder(1)
	Print FullDataFolder
	
	if (StrSearch(FullDataFolder, "In_Phase", 2) != -1)
		WaveEnding = "xx"
		Phase = "In Phase"
	elseif (StrSearch(FullDataFolder, "Out_of_Phase", 2) != -1)
		WaveEnding = "yx"
		Phase = "Out of Phase"
	elseif (StrSearch(FullDataFolder, "RightHand", 2) != -1)
		WaveEnding = "r"
		Phase = "Right Hand"
	elseif (StrSearch(FullDataFolder, "LeftHand", 2) != -1)
		WaveEnding = "l"
		Phase = "Left Hand"
	elseif (StrSearch(FullDataFolder, "Faraday", 2) != -1)
		WaveEnding = ""
		Phase = "Faraday Rotation"
	else
		abort
	endif
	
	String DataFolder
	DataFolder = GetDataFolder(0)
	
	DataFolder = ReplaceString(" ", DataFolder, "")
	DataFolder = ReplaceString("'", DataFolder, "")
	DataFolder = ReplaceString("-", DataFolder, "Minus")
	DataFolder = ReplaceString(".", DataFolder, "p")
	
	if(cmpstr("Mag_Transmissions", DataFolder) == 0)
		DataFolder = "MagTrans"
	endif
	
	if (StrSearch(FullDataFolder, "Constant_Field", 2) != -1)
		WaveTItle = "Field" + DataFolder
		AxisTItle = "Field " + DataFolder
		
	else
		AxisTitle = DataFolder + WaveEnding
		WaveTitle = DataFolder + WaveEnding
		Print WaveTitle
	endif
	
	//---Lets kill the window so we can recreate it, doesn't throw an error if the window doesn't exist yet---//
	DoWindow/K $(WaveTitle + "_vs_Frequency")
	//---Now let's create the window so we have make the graph---//
	Display/N= $(WaveTitle + "_vs_Frequency")
	
	String objName
		Variable index = 0
		do
			objName = GetIndexedObjName(":", 1, index)
			if (strlen(objName) == 0)
				break
			endif
			
			AppendToGraph/W =  $(WaveTitle + "_vs_Frequency") $ObjName
			ModifyGraph/Z/w=$(WaveTitle + "_vs_Frequency") offset($ObjName)={0, OSA*index}
			
			//---Make it look like we want it to---///
			ModifyGraph/Z cmplxMode=3, lsize=2, mirror=2, fStyle=1, fSize=14, axThick=2, nticks=10;
			Setaxis/Z/W= $('WaveTitle' + "_vs_Frequency") bottom 0.2, 2.25;
			Setaxis/Z/W= $('WaveTitle' + "_vs_Frequency")/A=2 left; 
			//Legend/C/N=text0/F=0/H={30,1,10}/A=MC
			Label bottom "\\Z16\\f01Frequency (THz)"
			Label left "\\Z16\\f01" + AxisTItle
			TextBox/C/N=text1/A=MC "\\f01\\Z12" + Phase
			TextBox/C/N=text1/X=-35.00/Y=41.00
			
			KBColorizeTraces#KBColorTablePopMenuProc(WaveTitle +"_vs_Frequency",0,"BlueRedGreen")	
			
			index += 1
		while(1)
		
End


function index_calc(tester_string)
	
	String tester_string
	Wave/C tester = $tester_string
	Variable delta = DimDelta(tester,0)
	//Print delta
	Variable i;
	//Variable d = 10
	
	Make/O/N=(numpnts(tester)), transmission
	transmission = cabs(tester);
	SetScale/P x 0,(delta),"", transmission
	//smooth 16, transmission
	
	Make/O/N=(numpnts(tester)), phase
	phase = Imag(r2polar(tester));
	Unwrap 2*pi, phase
	SetScale/P x 0,(delta),"", phase
	
	Make/O/N=(numpnts(tester)), freq
	for(i = 0; i<numpnts(freq); i+=1)
		freq[i] = delta*i
	endfor
	SetScale/P x 0,(delta),"", freq
end

function TDNR(trans, phi, f, d, num)
// This is a numerical procedure to solve (n,k) from the complex transmission of terahertz beam from thick single crystal with no substrate.
// This procedure defines a function TDNR, which takes the input of transmission amplitude, phase(unwrapped), frequency of the transmitted wave, 
// as well as the thickness of the crystal in units of 0.1mm, and number of iterations (usually 10-15 is more than enough, 20 would be overkill);
// the output is a series of waves, including n and k. Delta_T and Delta_Phi are the errors in the transmission amplitude and phase from the final values of n and k.

// This method is a 2D Newton method. It uses the so called "one at a time" stepper algorithm. The only simplification in the method is
// arctan(k/n) ~ k/n, which is valid except for small n and large k situation. Even then it only creates substantial error at low frequencies.  

// For this method to be valid, n*d need to be large enough so that the echoes from multiple reflections are well separated from the main transmission.

	wave trans, phi, f
	variable d, num
	wave/Z T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi, n, k
	variable i, c
	c=0.3
	
	if (f[0]==0)
	DeletePoints 0,1, trans,phi,f 
	endif
		
	make/o/n=(dimsize(phi,0)) n
	make/o/n=(dimsize(phi,0)) k
	make/o/n=(dimsize(phi,0)) T_n
	make/o/n=(dimsize(phi,0)) T_k
	make/o/n=(dimsize(phi,0)) Phi_n
	make/o/n=(dimsize(phi,0)) Phi_k
	make/o/n=(dimsize(phi,0)) Delta_T
	make/o/n=(dimsize(phi,0)) Delta_Phi
	
	n=phi*c/2/pi/d/f+1
	k=-c/2/pi/d/f*ln(trans*(n+1)^2/4/n)
	
	for(i=0;i<num;i+=1)
	
		Delta_T=(4*sqrt(n^2+k^2)/((n+1)^2+k^2)*exp(-k*2*pi*f*d/c)-trans)
		T_n=4*exp(-k*2*pi*f*d/c) / ((n+1)^2+k^2)^2 * (n*((n+1)^2+k^2)/sqrt(n^2+k^2)-sqrt(n^2+k^2)*2*(n+1))
		T_k=4*sqrt(n^2+k^2)/((n+1)^2+k^2)*exp(-k*2*pi*f*d/c)*(-2*pi*f*d/c)+4*exp(-k*2*pi*f*d/c)/((n+1)^2+k^2)^2*(k*((n+1)^2+k^2)/sqrt(n^2+k^2)-sqrt(n^2+k^2)*2*k)
		n=n-Delta_T*T_n/(T_n^2+T_k^2)
		k=k-Delta_T*T_k/(T_n^2+T_k^2)
		Delta_Phi=((n-1)*2*pi*f*d/c+atan(k/n)-2*atan(k/(n+1))-phi)
		Phi_n=2*pi*f*d/c-k/n^2+2*k/(n+1)^2
		Phi_k=1/n-2/(n+1)
		n=n-Delta_Phi*Phi_n/(Phi_n^2+Phi_k^2)
		k=k-Delta_Phi*Phi_k/(Phi_n^2+Phi_k^2)
		//Print "Iteration" + num2str(i)
		
	endfor
end

function TDNR_Chris_nc(trans, phi, f, d, num)
// This code is written for a thick single crystal that's mounted on top of a thick substate.  nc and kc are the real and imaginary parts of the substrates index of refraction.

// This method is a 2D Newton method. It uses the so called "one at a time" stepper algorithm. The only simplification in the method is
// arctan(k/n) ~ k/n, which is valid except for small n and large k situation. Even then it only creates substantial error at low frequencies.  

// For this method to be valid, n*d need to be large enough so that the echoes from multiple reflections are well separated from the main transmission.
	
	wave trans, phi, f
	variable d, num
	wave/Z T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi, n, k
	variable i, c
	c=0.3
	
	//-----Substrate Index of Refraction------//
	NVAR ns = root:Packages:NickLTHzExtras:ns
	
	if (f[0]==0)
	DeletePoints 0,1, trans,phi,f 
	endif
		
	make/o/n=(dimsize(phi,0)) n
	make/o/n=(dimsize(phi,0)) k
	make/o/n=(dimsize(phi,0)) T_n
	make/o/n=(dimsize(phi,0)) T_k
	make/o/n=(dimsize(phi,0)) Phi_n
	make/o/n=(dimsize(phi,0)) Phi_k
	make/o/n=(dimsize(phi,0)) Delta_T
	make/o/n=(dimsize(phi,0)) Delta_Phi
	make/o/n=(dimsize(phi,0)) w
	
	n=phi*c/2/pi/d/f+1
	k=-c/2/pi/d/f*ln(trans*(n+1)^2/4/n)
	w = 2*pi*f
	Variable nc=ns;
	Variable kc = 0;
	
	for(i=0;i<num;i+=1)
	
		Delta_T = (2*exp(-((d*k*w)/c))*Sqrt(k^2 + n^2)*Sqrt(kc^2 + (1 + nc)^2))/(Sqrt(k^2 + (1 + n)^2)*Sqrt((k + kc)^2 + (n + nc)^2)) - trans
		T_n = -(2*exp(-((d*k*w)/c))*Sqrt(kc^2 + (1 + nc)^2)*(2*k^3*kc - 2*k*kc*n*(1 + n) + k^4*(1 + n + nc) + k^2*(kc^2 + 2*n^2*(1 + n) + nc + 2*n*(2 + n)*nc + nc^2) + n*(1 + n)*(-kc^2 + (n^2 - nc)*(n + nc))))/(Sqrt(k^2 + n^2)*(k^2 + (1 + n)^2)^(3/2)*((k + kc)^2 + (n + nc)^2)^(3/2))
		T_k = -(2*exp(-((d*k*w)/c))*Sqrt(kc^2 + (1 + nc)^2)*(c*(k^5 + k^4*kc + 2*k^3*n^2 + kc*n^2*(1 + n)^2 + k^2*kc*(-1 + 2*(-1 + n)*n) - k*(-n^4 + kc^2*(1 + 2*n) + 4*n^2*nc + nc^2 + 2*n*nc*(1 + nc))) + d*(k^2 + n^2)*(k^2 + (1 + n)^2)*((k + kc)^2 + (n + nc)^2)*w))/(c*Sqrt(k^2 + n^2)*(k^2 + (1 + n)^2)^(3/2)*((k + kc)^2 + (n + nc)^2)^(3/2))
		n=n-Delta_T*T_n/(T_n^2+T_k^2)
		k=k-Delta_T*T_k/(T_n^2+T_k^2)
		Delta_Phi=(d*(-1 + n)*w)/c + atan(k/n) - atan(k/(1 + n)) + atan(kc/(1 + nc)) - atan((k + kc)/(n + nc))-phi
		Phi_n = -(k/(k^2 + n^2)) + k/(k^2 + (1 + n)^2) + (k + kc)/((k + kc)^2 + (n + nc)^2) + (d*w)/c
		Phi_k = n/(k^2 + n^2) - (1 + n)/(k^2 + (1 + n)^2) - (n + nc)/((k + kc)^2 + (n + nc)^2)
		n=n-Delta_Phi*Phi_n/(Phi_n^2+Phi_k^2)
		k=k-Delta_Phi*Phi_k/(Phi_n^2+Phi_k^2)
		
	endfor

end

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//	Below are the functions which control each of the buttons and pop up menus in the Analysis Palette.
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function FolderGrab(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			GetFolder()
			
		elseif(j==1 && i ==0)
			GetFolder()
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
			
		endif
		
	endswitch

	return 0
End

Function RefsInPhase(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			LoadScans("In_Phase", "Ref_Name")
			
		elseif(j==1 && i ==0)
			LoadScans("In_Phase", "Ref_Name")
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
			
		endif
	endswitch
	return 0
End

Function RefsOutofPhase(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				LoadScans("Out_Of_Phase", "Ref_Name")
			else
				Print "You must use the rotator to load out of phase scans"
			endif
		elseif(j==1 && i ==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				LoadScans("Out_Of_Phase", "Ref_Name")
			else
				Print "You must use the rotator to load out of phase scans"
			endif
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
			
		endif
			
	
	endswitch

	return 0
End

Function SampsInPhase(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			LoadScans("In_Phase", "Sample_Name")
			
		elseif(j==1 && i ==0)
			LoadScans("In_Phase", "Sample_Name")
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
			
		endif
	
	endswitch

	return 0
End

Function SampsOutPhase(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				LoadScans("Out_of_Phase", "Sample_Name")
			else
				Print "You must use the rotator to load out of phase scans"
			endif
			
		elseif(j==1 && i ==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				LoadScans("Out_of_Phase", "Sample_Name")
			else
				Print "You must use the rotator to load out of phase scans"
			endif
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
			
		endif
	
	endswitch

	return 0
End

Function TransCalc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				CalcTrans("In_Phase")
				CalcTrans("Out_of_Phase")
				MakePlot(2, "In_Phase", "Transmissions", "", "Txx")
				MakePlot(2, "Out_Of_Phase", "Transmissions", "", "Tyx")
			else
				CalcTrans("In_Phase")
				MakePlot(2, "In_Phase", "Transmissions", "", "Txx")
			endif
			
		elseif(j==1 && i ==0)
			ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
			If (V_Value==1)
				CalcTrans("In_Phase")
				CalcTrans("Out_of_Phase")

			else
				CalcTrans("In_Phase")
	
			endif
			
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
	
	endswitch

	return 0
End

Function ConvertToCircBasis(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			CircBasisTrans()
			
		elseif(j==1 && i ==0)
			CircBasisTrans()
	
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
			
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
	
	endswitch

	return 0
End

Function ConductivityCalc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
				If (PopVal == 1)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						SCCond()
						CircBasisSCCond()
					else
						SCCond()
						CircBasisSCCond()
						MakePlot(2, "In_Phase", "Sigma1", "", "Sigma1xx")
						MakePlot(2, "In_Phase", "Sigma2", "", "Sigma2xx")
					endif
				elseif (PopVal == 2)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						SCCond()
						CircBasisSCCond()
					else
						SCCond()
						CircBasisSCCond()
						MakePlot(2, "In_Phase", "Sigma1", "", "Sigma1xx")
						MakePlot(2, "In_Phase", "Sigma2", "", "Sigma2xx")
					endif
				elseif (PopVal == 3)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						TFCond()
						CircBasisTFCond()
					else
						TFCond()
						CircBasisTFCond()
						MakePlot(2, "In_Phase", "Conductance1", "", "G1xx")
						MakePlot(2, "In_Phase", "Conductance2", "", "G2xx")
					endif
				endif
		elseif(j==1 && i ==0)
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
				If (PopVal == 1)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						SCCond()
						CircBasisSCCond()
						
					else
						SCCond()
						CircBasisSCCond()
				
					endif
				elseif (PopVal == 2)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						SCCond()
						CircBasisSCCond()
						
					else
						SCCond()
						CircBasisSCCond()
	
					endif
				elseif (PopVal == 3)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)					
						TFCond()
						CircBasisTFCond()
					else
						TFCond()
						CircBasisTFCond()

					endif
				endif
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
		
	endswitch

	return 0
End

Function SuceptibilityCalc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
				If (PopVal == 1)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)
						Susceptibility()
						CircBasisSusceptibility()

					else
						Susceptibility()
						CircBasisSusceptibility()
						MakePlot(2, "In_Phase", "Chi1", "", "Chi1xx")
						MakePlot(2, "In_Phase", "Chi2", "", "Chi2xx")

					endif
					
				elseif (PopVal == 2)
					ControlInfo/W=NickL_THzAnalysis_Palette Rotator_CheckBox
					if (V_Value == 1)
						Susceptibility()
						CircBasisSusceptibility()

					else
						Susceptibility()
						CircBasisSusceptibility()
						MakePlot(2, "In_Phase", "Chi1", "", "Chi1xx")
						MakePlot(2, "In_Phase", "Chi2", "", "Chi2xx")

					endif
					
				elseif (PopVal == 3)
					Print "Susceptibility calculation can only be done on single crystal samples."
					
				endif
				
		elseif(j==1 && i ==0)
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
				If (PopVal == 1)
					Susceptibility()
					CircBasisSusceptibility()
					
				elseif (PopVal == 2)
					Susceptibility()
					CircBasisSusceptibility()
					
				elseif (PopVal == 3)
					Print "Susceptibility calculation can only be done on single crystal samples."
					
				endif
				
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
		
	endswitch

	return 0
End

Function FaradayCalc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			FaradayRotation()
		elseif(j==1 && i ==0)
			FaradayRotation()
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
		
	endswitch

	return 0
End

Function ConstFieldWaves(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			print "Constant field waves can only be generated on the magnet!"
		elseif(j==1 && i ==0)
			GenerateConstantField()
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
		
	endswitch

	return 0
End

Function FrequencyCuts(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up	
		
		variable i,j;
		ControlInfo/W=NickL_THzAnalysis_Palette BLC_CheckBox
		i = V_Value
		ControlInfo/W=NickL_THzAnalysis_Palette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
				If (PopVal == 1)
					FreqCutsSCZeroField()
				elseif (PopVal == 2)
					FreqCutsSCZeroField()
				elseif (PopVal == 3)
					FreqCutsSCZeroField()
				endif
		elseif(j==1 && i ==0)
			//CircBasis()
		elseif(i==1 && j==1)
			print "You must select either the BLC or the Magnet, not both!"
		elseif(i==0 && j==0)
			print "You must select either the BLC or the Magnet!"
		endif
	endswitch

	return 0
End

Function MeasurementType(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			SVAR SampType = root:Packages:NickLTHzExtras:SampType
			NVAR PopVal = root:Packages:NickLTHzExtras:PopVal
			
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			PopVal = PopNum
			SampType = popStr
			break
	endswitch

	return 0
End

Function OffsetType(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
		
			SVAR OsetType = root:Packages:NickLTHzExtras:OsetType
			
			String PopStr = pa.popStr
			OsetType = PopStr
			
			break
	endswitch

	return 0
End

Function MakeAGraph(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			
			GenerateGraph()
	
	endswitch

	return 0
End


Function CheckGraphs(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Rotator(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Magnet(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BLC(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Subn(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

