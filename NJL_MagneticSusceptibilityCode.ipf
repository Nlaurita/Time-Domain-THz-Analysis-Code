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
//  This is a procedure, written by Dr. Nicholas J. Laurita on July 13th, 2017, which calculates the magnetic susceptibility of a single crystal magnetic insulator from the complex transmission spectrum
//  obtained by time-domain terahertz spectroscopy.  This is the most general form of this code and should be easily adaptable to any group who wishes to calculate the magnetic susceptiblity of a sample
//  from TDTS data.
//
//  Much more information on the analysis routine for computing the magnetic susceptibility in a time-domain THz measurement can be found in my Ph.D. thesis which can be found here:
//  https://drive.google.com/file/d/0B7K_8wnuzg_QRE9lQkczZlpEUDQ/view
//
//*******************************************************************************************************************************************************************************************************************************************************************//



//************************************************************************************************** INSTRUCTIONS FOR RUNNING THE CODE *************************************************************************************************************//
//
// 		Press and hold Control-2 to pull up the magnetic susceptiblity analysis palette.
//
// #1	Select the type of experiment you have performed.  The code can analyze data taken as a function of temperature or magnetic field, although the choice only changes the naming convention the program looks 
//		for and could easily be adapted to analyze data taken as a function of any parameter.  
//
// #2	Enter in the thickness of your sample in millimeters into the "Sample Thickness" text box.
//
// #3	Enter your sample name in the "Sample Name" text box on the panel.  All the data you wish to analyze must be in a folder in the data browser named "Sample_Name"_Analysis.  For instance,
//		If you were analyzing data taken on FeSc2S4, then you would write "FeSc2S4" in the text box in the panel and then all data must be stored in a folder called "FeSc2S4_Analysis" for the program to run.
//
// #4	The appropriate text waves and folders must be created before calculating the magnetic susceptibility.  These textwaves tell the code which folders to look into for the transmission and which folders to create 
//		during the analysis. All of the following waves and folders must be contained in the master "SampleName_Analysis" folder.
//
//				For temperature dependent data:
//					Create a TEXT wave called "Temps" which contains the temperatures at which data was taken.  For instance elements in this wave might be "4K", "10K", etc.  Create a folder called 
//					"transmissions" where your complex transmission waves will be.  In this folder you must name each transmission in the form: "T_SampleName_temperature".  For example,
//					if you have measured FeSc2S4 at 10K and 20K then you would have complex transmission waves in the "transmissions" folder with the names "T_FeSc2S4_10K" and 
//					"T_FeSc2S4_20K".  You would also have elements in the text wave "Temps" of "10K" and "20K".
//
//				For magnetic field dependent data:
//					Create a TEXT wave called "Temps" which contains the temperatures at which taken was taken.  For instance elements might be "4K", "10K", etc.  Additionally you must create a TEXT wave
//					for each temperature with the name "Fields_Temps" that contains the fields that were measured at each temperature.  For instance if data was taken at 15K then there will be an additional text 
//					wave with the name "Fields_15K".  Populate this text wave with the magnetic fields which were measured at 15K.  The format requires two digits after the decimal point, so elements 
//					could be"1.50" or "15.00" etc. No units are neeeded here but the wave must still be a text wave.
//
//					Inside the "SampleName_Analysis" folder, create a new folder for each temperature.  Again these folders will have names like "4K", "10K" etc.  Then, in each of these folders then create a new folder named
//					"transmissions" and load your complex transmission waves into these folders.  The naming convention for magnetic field dependent data is "T_Temps_FieldskG".  For example, if you have measured
//					FeSc2S4 at 10K at fields of 20 kG and 30 kG (We measure our magnetic fields in killoGauss) then you would have complex transmission waves with the titles "T_10K_20.00kG" and "T_10K_30.00kG"
//					in your "transmissions" folder, which is in your "10K" folder, which is in your "FeSc2S4_Analysis" folder.  Do this for every temperature and field. 
//
// #5 	Once your data is loaded in the correct format and you've created your temperature and field text waves then you're almost ready to calclulate the magnetic susceptibility.  The last thing to do is select a 
//		reference temperature and enter it into the "Reference Temp:" text box.  An ideal reference temperature is a temperature at which the dielectric contributions to the index of refraction will not appreciably change below
//		 AND magnetic correlations have not yet onset.  I would suggest a temperature above a magnetic transition or above the appearance of a magnetic excitation but below roughly 100K.  For instance, you may chose 80K
//		 as a reference temperature and then should inpurt "80K" into the text box.  For more information regarding why a reference temperature is needed please see Ch. 2 of my Ph.D. thesis.  
//
// 		Once you've selected a reference temperature, hit the "Calculate Magnetic Susceptibilty" button.  The code will take your transmission, and calculate the real and imaginary parts of the index of refraction, 
//		the phase, and the real and imagainary parts of the complex magnetic susceptibility.
//
// #6	I've also included the "Make Graph" feature in this analysis palette because it's a very useful function.  Set any data folder as the default folder in the data brower and then hit the "Make Graph"
//		button to automatically generate a publishing ready plot of the data in that folder.  Its a great way to plot a lot of data very quickly.  The "Offset Amount" option creates an offset between waves in the plot to 
//		better display the data.
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
	DoWindow/HIDE=? $("SusceptbilityPalette")
	if (V_flag != 0)
		DoWindow/F SusceptbilityPalette;
	else
		Execute/Q "SusceptbilityPalette()"
	endif
	
	//---Lets define some global variables---//

	//---Folder where data will be loaded from.  It gets generated by the "Got Folder Info" button in the analysis palette---//
	String/G base_folder
	
	//---Sample_Name stores the name to look for in the sample scans.  For instance in "SmB6_5K_1", Sample_Name = SmB6---//
	String/G Sample_Name = strVarOrDefault("root:Packages:NickLTHzExtras:Sample_Name", "nan") 
	
	//---dsc is the single crystal sample thickness, thickness has to be in mm's.---//
	variable/G dsc = NumVarOrDefault("root:Packages:NickLTHzExtras:dsc", 0) 
	
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
Window SusceptbilityPalette() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1494,316,1876,703)
	ModifyPanel frameStyle=3
	SetDrawLayer UserBack
	SetDrawEnv fsize= 16,fstyle= 1,textyjust= 2
	DrawText 32,3,"Magnetic Susceptibility Analysis Palette"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,334,375,334
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 31,358,"Please Read Instruction In Code Before Use"
	SetDrawEnv fsize= 14
	DrawText 21,374,"You can use CRTL-2 to pull up the analysis palette"
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 66,40,"Data Was Taken As A Function Of:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,24,375,24
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,133,375,133
	SetDrawEnv fsize= 14,fstyle= 1,textyjust= 1
	DrawText 118,147,"Enter Sample Name"
	SetDrawEnv fsize= 14,fstyle= 1,textxjust= 2,textyjust= 1
	DrawText 271,93,"Enter Sample Thicnkess"
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 97,280,"Make Graph / Add Offset:"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,80,375,80
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 74,206,"Calculate Magnetic Susceptibility"
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,186,375,186
	SetDrawEnv linethick= 3,linebgc= (56576,56576,56576)
	DrawLine 0,255,375,255
	SetVariable SampleName,pos={5,159},size={204,20},title="\\Z14\\F'Arial'Sample Name:     "
	SetVariable SampleName,value= root:Packages:NickLTHzExtras:Sample_Name
	SetVariable SampleThickness,pos={5,105},size={200,20},title="\\Z14\\F'Arial'Sample Thickness (mm)"
	SetVariable SampleThickness,value= root:Packages:NickLTHzExtras:dsc
	SetVariable OffsetAmount,pos={215,296},size={153,20},title="\\Z14\\F'Arial'Offset Amount"
	SetVariable OffsetAmount,value= root:Packages:NickLTHzExtras:OSA
	CheckBox BLC_Checkbox,pos={35,56},size={100,16},proc=BLC,title="\\Z14\\F'Arial'Temperature"
	CheckBox BLC_Checkbox,value= 0
	CheckBox DoGraphs_Checkbox2,pos={-69,-100},size={129,16},proc=CheckGraphs,title="\\Z14\\F'Arial'Generate Graphs"
	CheckBox DoGraphs_Checkbox2,value= 1
	CheckBox Magnet_Checkbox,pos={227,56},size={110,16},proc=Magnet,title="\\Z14\\F'Arial'Magnetic Field"
	CheckBox Magnet_Checkbox,value= 0
	Button button7,pos={4,210},size={180,40},proc=SuceptibilityCalc,title="Calculate Magnetic \rSusceptibility"
	Button button7,fSize=12,fStyle=1
	SetVariable RefTemp,pos={204,219},size={157,20},title="\\F'Arial'\\Z14Reference Temp"
	SetVariable RefTemp,value= root:Packages:NickLTHzExtras:RefTemp
	Button button0,pos={3,287},size={180,40},proc=MakeAGraph,title="Make Graph"
	Button button0,fSize=12,fStyle=1
EndMacro

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function calculates the Susceptibility for a single crystal for either the BLC or Magnet, with or without the rotator.  Grabs an index of refraction at some higher temperature above the magnetism in the sample which
//  the user specifies for the calculation.  This temperature is call Ref_Temp.  Note there is no function for calculating the susceptibility of a thin film.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
Function Susceptibility()

//-----String where you enter in the name of the sample-----//
	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name

//-----Thickness of the sample in mm's------//
	NVAR dsc = root:Packages:NickLTHzExtras:dsc

//-----Reference Temperature for Susceptibility Calculation------//
	SVAR RefTemp = root:Packages:NickLTHzExtras:RefTemp
	
//------Declare other variables and dummy waves------//
	variable delta, i, j, c = 299792458, mu = 1//(1.25663706*10^-6);
	Variable num_iter = 10  //--The number of Netwon-Rhapson iterations
	
//-----Dummy variables-----//
	String fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8, fname9
	
	//---This is the case that we're analyzing data as a function of temperature---//
	ControlInfo/W=SusceptbilityPalette BLC_CheckBox
		if (V_Value == 1)
		
			//---Sets the data folder to the Analysis folder just in case---//
			SetDataFolder root:$Sample_Name + "_Analysis"
			
			//---recognizes the wave that we'll need while running this function---//
			Wave/T Temps
			 Variable NumTemps = numpnts(Temps)
			
					//---Creates folders for the complex susceptibility, complex n, and phase---//
					NewDataFolder/O root:$Sample_Name + "_Analysis":n
					NewDataFolder/O root:$Sample_Name + "_Analysis":k
					NewDataFolder/O root:$Sample_Name + "_Analysis":Phase
					NewDataFolder/O root:$Sample_Name + "_Analysis":Chi1
					NewDataFolder/O root:$Sample_Name + "_Analysis":Chi2
					
					for(i = 0; i < NumTemps; i+=1)
						//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
						SetDataFolder root:$Sample_Name + "_Analysis":Transmissions
						index_calc("T_" + Sample_Name + "_" + Temps[i])
						Wave transmission, phase, freq
						
						//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
						CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
						Wave W_coef
						phase -= W_coef[0]
						
						//-----Run the newton-raphson code to find the index-----//
						delta = DimDelta($("T_" + Sample_Name + "_" + Temps[i]),0)
						
						//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
						TDNR(transmission, phase, freq, dsc, num_iter)
						
						//---Declare n and k and then set their scale---//
						Wave n,k
						SetScale/P x 0,(delta),"", n,k
						
						//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
						Duplicate/O $("n"), root:$Sample_Name + "_Analysis":n:$("n_" + Sample_Name + "_" + Temps[i])
						KillWaves $("n")
						
						Duplicate/O $("k"), root:$Sample_Name + "_Analysis":n:$("k_" + Sample_Name + "_" + Temps[i])
						Duplicate/O $("k"), root:$Sample_Name + "_Analysis":k:$("k_" + Sample_Name + "_" + Temps[i])
						KillWaves $("k")
						
						Duplicate/O $("phase"), root:$Sample_Name + "_Analysis":Phase:$("Theta_" + Sample_Name + "_" + Temps[i])
						KillWaves $("phase")
						
					endfor
					
						//-----A bit of wave clean up-----//
						Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
						KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					
					for(i = 0; i < numpnts(Temps); i+=1)
						
						//---Now go to the Mag Transmissions folder so we can calculate the complex susceptibility.  n, k, phase, and MagT should all be there already---//
						SetDataFolder root:$Sample_Name + "_Analysis":n
						
						//---Name some strings for waves we'll need later---//
						fname1 = ("n_" + Sample_Name + "_" + Temps[i])
						fname2 = ("k_" + Sample_Name + "_" + Temps[i])
						
						fname3 = ("n_" + Sample_Name + "_" + RefTemp)
						fname4 = ("k_" + Sample_Name + "_" + RefTemp)
						
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname6 = ("nComplex_" + Sample_Name + "_" + RefTemp)
						
						fname7 = ("Chi2_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1_" + Sample_Name + "_" + Temps[i])
						
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
						
						F6 = F3 + sqrt(-1)*F4
						F5 = F1 + sqrt(-1)*F2
						
						F9 = (F5/F6)^2 - 1
						
						F8 = real(F9)
						F7 = imag(F9)
						
					endfor
					
					//---A loop which will do some wave clean up for us---//
					for(i = 0; i < numpnts(Temps); i+=1)
						
						fname2 = ("k_" + Sample_Name + "_" + Temps[i])
						fname5 = ("nComplex_" + Sample_Name + "_" + Temps[i])
						fname7 = ("Chi2_" + Sample_Name + "_" + Temps[i])
						fname8 = ("Chi1_" + Sample_Name + "_" + Temps[i])
						fname9 = ("ComplexChi_" + Sample_Name + "_" + Temps[i])
						
						Duplicate/O $fname7, root:$Sample_Name + "_Analysis":Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_Analysis":Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
					endfor
		endif
		
	//---This is the case that we're running the Magnet---//
	ControlInfo/W=SusceptbilityPalette Magnet_CheckBox
	if (V_Value == 1)
	
		SetDataFolder root:$Sample_Name + "_Analysis":
		Wave/T Temps
		NumTemps = NumPnts(Temps)
		
		for (i=0; i<NumTemps; i+=1)
		
			//---Set the data folder to the main folder to find "Fields_Temps" waves
			SetDataFolder root:$Sample_Name + "_Analysis"
		
			//---Creates folders for the complex conductivity, complex n, and phase---//
			NewDataFolder/O root:$Sample_Name + "_Analysis":$Temps[i]:n
			NewDataFolder/O root:$Sample_Name + "_Analysis":$Temps[i]:k
			NewDataFolder/O root:$Sample_Name + "_Analysis":$Temps[i]:Phase
			NewDataFolder/O root:$Sample_Name + "_Analysis":$Temps[i]:Chi1
			NewDataFolder/O root:$Sample_Name + "_Analysis":$Temps[i]:Chi2
			
			string fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			Variable NumFields = NumPnts($fname)
		
			//--- j runs over the number of Fields we have---//
			for (j=0; j<NumFields; j+=1)
		
				//----Create real versions of all the waves, the Newton Raphson needs non-complex waves----//
				SetDataFolder root:$Sample_Name + "_Analysis":$Temps[i]:Transmissions
				index_calc("T_" + Temps[i] + "_" + Fields[j] + "kG")
				Wave transmission, phase, freq
				
				//------Corrects if there is an offset in the phase, should intersect at 0, not 2 pi------//
				CurveFit/Q/NTHR=0 line  phase(0.2,0.6)
				Wave W_coef
				phase -= W_coef[0]
				
				//-----Run the newton-raphson code to find the index-----//
				delta = DimDelta($("T_" + Temps[i] + "_" + Fields[j] + "kG"),0)
				
				//---Normal Newton_Raphson code for a single crystal with an aperture reference---//
				TDNR(transmission, phase, freq, dsc, num_iter)
				
				//---Declare n and k and then set their scale---//
				Wave n,k
				SetScale/P x 0,(delta),"", n,k
				
				//-----Sort n, k, phase into the appropriate folders labeled by the current temperature-----//
				Duplicate/O $("n"), root:$Sample_Name + "_Analysis":$Temps[i]:n:$("n_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("n")
				
				Duplicate/O $("k"), root:$Sample_Name + "_Analysis":$Temps[i]:n:$("k_" + Temps[i] + "_" + Fields[j] + "kG")
				Duplicate/O $("k"), root:$Sample_Name + "_Analysis":$Temps[i]:k:$("k_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("k")
				
				Duplicate/O $("phase"), root:$Sample_Name + "_Analysis":$Temps[i]:Phase:$("Theta_" + Temps[i] + "_" + Fields[j] + "kG")
				KillWaves $("phase")
					
			endfor
				
					//-----A bit of wave clean up-----//
					Wave transmission, phase, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
					KillWaves/Z transmission, phase, freq, W_coef, W_sigma, n, k, T_n, T_k, Phi_n, Phi_k, Delta_T, Delta_Phi
		endfor
		
		for (i=0; i<NumTemps; i+=1)	
			SetDataFolder root:$Sample_Name + "_Analysis"
					
			fname = "Fields_" + Temps[i]
			Wave/T Fields = $fname
			NumFields = NumPnts($fname)
			
			for (j=0; j<NumFields; j+=1)
				
				//---Name some strings for waves we'll need later---//
				fname1 = ("n_" + Temps[i] + "_" + Fields[j] + "kG")
				fname2 = ("k_" + Temps[i] + "_" + Fields[j] + "kG")
				fname3 = ("n_" + RefTemp + "_" + Fields[j] + "kG")
				fname4 = ("k_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
				fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
				
				fname7 = ("Chi2_" + Temps[i] + "_" + Fields[j] + "kG")
				fname8 = ("Chi1_" +  Temps[i] + "_" + Fields[j] + "kG")
				
				fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
				
				SetDataFolder root:$Sample_Name + "_Analysis":$(RefTemp):n
				Duplicate/O $fname3, root:$Sample_Name + "_Analysis":$Temps[i]:n:$fname3
				Duplicate/O $fname4, root:$Sample_Name + "_Analysis":$Temps[i]:n:$fname4
				
				SetDataFolder root:$Sample_Name + "_Analysis":$Temps[i]:n
				
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
				
				F6 = F3 + sqrt(-1)*F4
				F5 = F1 + sqrt(-1)*F2
				
				F9 = (F5/F6)^2 - 1
				
				F8 = real(F9)
				F7 = imag(F9)
		
			endfor
				
				//---A loop which will do some wave clean up for us---//
				
				if (cmpstr(Temps[i], RefTemp) == 0)
					for (j=0; j<NumFields; j+=1)
					
						fname2 = ("k_" + Temps[i] + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname7 = ("Chi2_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_Analysis":$Temps[i]:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_Analysis":$Temps[i]:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname5, $fname7, $fname8, $fname9
						
					endfor
					
				else
					for (j=0; j<NumFields; j+=1)
						fname2 = ("k_" + Temps[i] + "_" + Fields[j] + "kG")
						fname3 = ("n_" + RefTemp + "_" + Fields[j] + "kG")
						fname4 = ("k_" + RefTemp + "_" + Fields[j] + "kG")
						fname5 = ("nComplex_" + Temps[i] + "_" + Fields[j] + "kG")
						fname6 = ("nComplex_" + RefTemp + "_" + Fields[j] + "kG")
						fname7 = ("Chi2_" + Temps[i] + "_" + Fields[j] + "kG")
						fname8 = ("Chi1_" +  Temps[i] + "_" + Fields[j] + "kG")
						fname9 = ("ComplexChi_" + Temps[i] + "_" + Fields[j] + "kG")
						
						Duplicate/O $fname7, root:$Sample_Name + "_Analysis":$Temps[i]:Chi2:$fname7
						Duplicate/O $fname8, root:$Sample_Name + "_Analysis":$Temps[i]:Chi1:$fname8
						
						Killwaves/Z $fname2, $fname3, $fname4, $fname5, $fname6, $fname7, $fname8, $fname9
					endfor
				endif
			endfor
	endif
End


////-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//// This function takes in as parameters the folder and the general wave name and then populates a graph of the data.  Should work for anything that's calculated via this panel in both the BLC or Magnet with or without 
//// the rotator.  Q is a parameter which permutates the possible paths for the data folder.  Folder1 refers to the "" or "Out_Of_Phase" part of the wave path.  Folder2 is the next folder down.  
////You can trick the code into running without the rotator by putting the second folder in for folder1 and  then having nothing in folder2.  For instance MakePlot("Transmissions", "", "T") would run without the rotator.
//// I added an additional folder, folder3, so as to make the code more versitile.  
////-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//
//Function MakePlot(Q, Folder1, Folder2, Folder3, WaveTitle)
//	//---Defines the number of folders we want to have in the title---//
//	Variable Q;
//	
//	//---Define the folder and the wave that we want on each graph---//
//	String WaveTItle, Folder1, Folder2, Folder3;
//	
//	//-----Folder where data is loaded from-----//
//	SVAR base_folder = root:Packages:NickLTHzExtras:base_folder
//	
//	//-----String where you enter in the name of the sample and aperature-----//
//	SVAR Sample_Name = root:Packages:NickLTHzExtras:Sample_Name
//	SVAR Ref_Name = root:Packages:NickLTHzExtras:Ref_Name
//	
//	//---Define some dummy variables---//
//	String fname, fname2, fname3, fname4;
//	variable i,j,k;
//	
//	//---Only create the graphs if the check box on the front panel is selected---//
//		ControlInfo/W=SusceptbilityPalette DoGraphs_CheckBox
//			if (V_Value == 1)
//	
//			//---This is the case that the "BLC" option is selected on the front panel---//
//				ControlInfo/W=SusceptbilityPalette BLC_CheckBox
//					if (V_Value == 1)
//	
//						//---recognizes the waves that we'll need while running this function---//
//						SetDataFolder root:$Sample_Name + "_Analysis"
//						Wave/T Temps
//						Wave num_samp_scans
//						Variable NumTemps = numpnts(Temps)
//						
//						if (Q==1)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1)
//						elseif (Q==2)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1):$(Folder2)
//						elseif (Q==3)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1):$(Folder2):$(Folder3)
//						else
//							abort
//						endif
//							
//							//---Lets kill the window so we can recreate it, doesn't throw an error if the window doesn't exist yet---//
//							DoWindow/K $(WaveTitle + "_vs_Frequency")
//							//---Now let's create the window so we have make the graph---//
//							Display/N= $(WaveTitle + "_vs_Frequency")
//							
//							//---Append the traces to the graph--//
//							for(i=0; i<numpnts(Temps); i+=1)
//								fname = (WaveTitle + "_" +  Sample_Name + "_" + Temps[i])
//								AppendToGraph/W =  $(WaveTitle + "_vs_Frequency") $fname
//							endfor
//							
//							//---Make it look like we want it to---///
//							ModifyGraph/Z cmplxMode=3, lsize=2, mirror=2, fStyle=1, fSize=14, axThick=2, nticks=10;
//							Setaxis/Z/W= $(WaveTitle + "_vs_Frequency") bottom 0.2, 2.25;
//							Setaxis/Z/W= $(WaveTitle + "_vs_Frequency")/A=2 left; 
//							Legend/C/N=text0/F=0/H={30,1,10}/A=MC
//							Label bottom "\\Z16\\f01Frequency (THz)"
//							Label left "\\Z16\\f01" + Folder1 + " " + Folder2
//								
//							KBColorizeTraces#KBColorTablePopMenuProc(WaveTitle +"_vs_Frequency",0,"BlueRedGreen")	
//					endif
//		
//		//---This is the case that we are analyzing data as a function of magnetic field--//
//		ControlInfo/W=SusceptbilityPalette Magnet_CheckBox
//			if (V_Value == 1)
//			
//				variable p=0, Num_TempFields;
//				
//				//--- This loop tells us how many temperatures were measured at in the Magnet system---//
//				for (k=1; k<=100; k+=1)
//					SetDataFolder root:$Sample_Name + "_Analysis":
//					fname = "TempFields_" + num2str(k)
//					Wave/T TFNow = $fname
//					p += WaveExists(TFNow)
//				endfor
//					Num_TempFields = p
//					
//					for (k=1; k<=Num_TempFields; k+=1)
//
//						SetDataFolder root:$Sample_Name + "_Analysis":
//					
//						String CurrentTF = "TempFields_" + num2str(k)
//						Wave/T TFNow = $CurrentTF
//						
//						Variable NumFields = numpnts(TFNow) 
//						
//						if (Q==1)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$TFNow[0]:$(Folder1)
//						elseif (Q==2)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1):$TFNow[0]:$(Folder2)
//						elseif (Q==3)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1):$TFNow[0]:$(Folder2):$(Folder3)
//						elseif (Q==4)
//							//---Sets the current folder---//
//							SetDataFolder root:$Sample_Name + "_Analysis":$(Folder1):$(Folder2):$TFNow[0]:$(Folder3)
//						else
//							abort
//						endif
//						
//						//---Lets kill the window so we can recreate it, doesn't throw an error if the window doesn't exist yet---//
//						DoWindow/K $(WaveTitle + "_vs_Frequency_" + TFNow[0])
//						//---Now let's create the window so we have make the graph---//
//						Display/N=$(WaveTitle + "_vs_Frequency_" + TFNow[0])
//
//						//---Loop that runs over fields---//
//						for( j = 1; j < NumFields; j += 1)
//
//							//---Append the traces to the graph--//
//								fname = (WaveTitle + "_" + TFNow[0] + "_" + TFNow[j] + "kG")
//								AppendToGraph/W = $(WaveTitle + "_vs_Frequency_" + TFNow[0]) $fname					
//							
//								//---Make it look like we want it to---///
//								ModifyGraph/Z cmplxMode=3, lsize=2, mirror=2, fStyle=1, fSize=14, axThick=2, nticks=10;
//								Setaxis/Z/W= $(WaveTitle + "_vs_Frequency_" + TFNow[0]) bottom 0.2, 2.25;
//								Setaxis/Z/W= $(WaveTitle + "_vs_Frequency_" + TFNow[0])/A=2 left; 
//								Legend/C/N=text0/F=0/H={30,1,10}/A=MC
//								Label bottom "\\Z16\\f01Frequency (THz)"
//								Label left "\\Z16\\f01" + Folder1 + " " + Folder2
//						
//								//---Colorize the trances so they look very nice---//
//								KBColorizeTraces#KBColorTablePopMenuProc(WaveTitle + "_vs_Frequency_" + TFNow[0],0,"BlueRedGreen")
//						endfor
//					endfor	
//			endif
//	endif
//End

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
// This function makes a plot of all the data that's contained in a given folder.  Simply select that folder as the default folder in the data browser and it will upload all the data to a plot.  Can add an offset as well. 
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function GenerateGraph()

	//-----OSA is the offset amount-----//
	NVAR OSA = root:Packages:NickLTHzExtras:OSA;
	
	//---Define some dummy variables---//
	String fname, fname2, fname3, fname4, WaveTitle, AxisTitle, WaveEnding, Phase
	variable i,j,k;
	
	String FullDataFolder
	FullDataFolder = GetDataFolder(1)
	
	String DataFolder
	DataFolder = GetDataFolder(0)
	
	DataFolder = ReplaceString(" ", DataFolder, "")
	DataFolder = ReplaceString("'", DataFolder, "")
	DataFolder = ReplaceString("-", DataFolder, "Minus")
	DataFolder = ReplaceString(".", DataFolder, "p")
	
	AxisTitle = DataFolder
	WaveTitle = DataFolder

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

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//	Below are the functions which control each of the buttons and pop up menus in the Analysis Palette.
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


Function SuceptibilityCalc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		
		variable i,j;
		ControlInfo/W=SusceptbilityPalette BLC_CheckBox
		i = V_Value
		ControlInfo/W=SusceptbilityPalette Magnet_CheckBox
		j = V_Value
		
		if(i==1 && j==0)
			Susceptibility()
		elseif(j==1 && i ==0)
			Susceptibility()	
		elseif(i==1 && j==1)
			print "You must select either temperature depenence or magnetic field dependence, not both."
		elseif(i==0 && j==0)
			print "You must select either temperature depenence or magnetic field dependence."
		endif
		
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



