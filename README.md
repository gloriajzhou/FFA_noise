**README**

Single-cell analysis reveals condition-dependent bioproduction noise and guides gene circuit design for reduced metabolite variations

The code runs using MATLAB software and has been tested on versions 2023b (https://www.mathworks.com/products/new_products/release2023b.html) and 2024b (https://www.mathworks.com/products/new_products/release2024b.html). 

**Instructions for Software Installation**

If your local device does not already have MATLAB installed, you can follow the links above for specific versions of the software, or search "download MATLAB" using Google. A MathWorks account and a valid license (or 30-day trial) will be required to download the installer. Specific system requirements can be found at https://www.mathworks.com/support/requirements/matlab-system-requirements.html. 

The code does not require non-standard hardware, and typical installation times can range from a few seconds to a few minutes. 

**Instructions for Running ngl_solv.m and nml_solv.m**

Download the MATLAB files (ngl_solv.m or nml_solv.m) to your local device. The code is self-contained and does not depend on any external datasets to run. Simulations can be tuned by adjusting the variables defined in the code. Confirm that the current working folder is set to the folder containing the code file that is being run (if it's not, MATLAB should prompt you to change the working folder after attempting to run the code).

The code should take a few seconds to run (less than a minute). After running the code, figures will automatically be generated showing relative product noise or relative mean product number as a function of promoter and repression strength (see Fig. 4c and 4d in the main text as examples). These figures can be saved to your local device by going to File > Save As. Datasets can be exported from the MATLAB environment by opening the Workspace pane and double-clicking the variable of interest (e.g., promoter_strength), which opens a new tab or window showing the specific data values. These values can then be copied and pasted into a spreadsheet (e.g., in Microsoft Excel) for further data processing. 

**Reproducing Manuscript Figures**

Directly downloading the code files and running them in MATLAB as-is will allow the quantitative data values in the manuscript to be reproduced (see Source Data file for Fig. 4c and 4d). See Supplementary Table 3 and Supplementary Table 4 in the Supplementary Information file for the references used to define the model parameters. Fig. 4c (NGL) and 4d (NML) in the main text were created by including both relative product noise and relative mean product number values in the same plot using primary and secondary axes in Microsoft Excel (version 16.107.3).
