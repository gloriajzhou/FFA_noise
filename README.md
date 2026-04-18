**README**

Single-cell analysis reveals condition-dependent bioproduction noise and guides gene circuit design for reduced metabolite variations

The code runs using MATLAB software and has been tested on versions 2023b (https://www.mathworks.com/products/new_products/release2023b.html) and 2024b (https://www.mathworks.com/products/new_products/release2024b.html). 

**Instructions for Software Installation**

If your local device does not already have MATLAB installed, you can follow the links above for specific versions of the software, or search "download MATLAB" using Google. A MathWorks account and a valid license (or 30-day trial) will be required to download the installer. Specific system requirements can be found at https://www.mathworks.com/support/requirements/matlab-system-requirements.html. 

The code does not require non-standard hardware, and typical installation times can range from a few seconds to a few minutes. 

**Instructions for Running ngl_solv.m and nml_solv.m**

Download the MATLAB files (ngl_solv.m or nml_solv.m) to your local device. The code is self-contained and does not depend on any external datasets to run. Simulations can be tuned by adjusting the variables defined in the code. Confirm that the current working folder is set to the folder containing the code file that is being run (if it's not, MATLAB should prompt you to change the working folder after attempting to run the code).

The code should take a few seconds to run (less than a minute). After running the code, figures will automatically be generated showing relative product noise or relative mean product number as a function of promoter and repression strength (see Fig. 4c and 4d in the main text as examples). These figures can be saved to your local device by going to File > Save As. Datasets can be exported from the MATLAB environment by opening the Workspace pane and double-clicking the variable of interest (e.g., promoter_strength), which opens a new tab or window showing the specific data values. These values can then be copied and pasted into a spreadsheet (e.g., in Microsoft Excel) for further data processing. 

**Instructions for Running eFAST Simulations**

Code files were adapted from previously published methods for the extended Fourier amplitude sensitivity test (eFAST) (see Hartline, C. J. & Zhang, F. The Growth Dependent Design Constraints of Transcription-Factor-Based Metabolite Biosensors. ACS Synth. Biol. 11, 2247–2258 (2022) and Marino, S., Hogue, I. B., Ray, C. J. & Kirschner, D. E. A Methodology For Performing Global Uncertainty And Sensitivity Analysis In Systems Biology. Journal of theoretical biology 254, 178 (2008)). 

Download the eFAST files folder to your local device. The code files within this folder are self-contained and do not depend on external datasets to run. The Parameter_settings.m files for NGL and NML specify model parameters, including the number of time points, the number of outputs, the number of parameters, the range of variations for each parameter, and the baseline values. efast_sd.m is a function file called in mainModel.m to calculate Si, Sti, rangeSi, and rangeSti values. modelODE.m for NGL and NML are files specifying the equations used to model each circuit. efast_ttest.m is a function file called in mainModel.m to calculate the p-values used for determining statistical significance. To run the sensitivity analysis, open mainModel.m for NGL or NML and click "Run" after confirming the desired parameter values and model specifications. The code should take less than a minute to run. The data for meanSti and stdSti can be found in the Workspace pane by double-clicking the parameters, while the p-values (p_Sti_out) can be obtained from the Command Window output. These values can be copied and pasted into a spreadsheet (e.g., in Microsoft Excel) for further data processing.

**Reproducing Manuscript Figures**

Directly downloading the code files and running them in MATLAB as-is will allow the quantitative data values in the manuscript to be reproduced for Fig. 4c and 4d in the Source Data file. Due to the stochastic nature of the eFAST model, the same command can generate new data with each simulation run from random sampling across the parameter space, so the exact values shown in Supplementary Fig. 6a-d can vary with re-runs of the same code, but the conclusions from the modeling data (e.g., relative mean values and statistical significance of each variable) remain unchanged. See Supplementary Table 3 and Supplementary Table 4 in the Supplementary Information file for the references used to define the model parameters in ngl_solv.m and nml_solv.m. Fig. 4c (NGL) and 4d (NML) in the main text were created by including both relative product noise and relative mean product number values in the same plot using primary and secondary axes in Microsoft Excel (version 16.107.3). Supplementary Fig. 6a-d were created by extracting the meanSti, stdSti, and p_Sti_out values for var1 (rel_mp or relative mean FFA) and var2 (rel_noise or relative FFA noise). 
