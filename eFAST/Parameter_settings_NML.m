%% PARAMETER INITIALIZATION

nTimePoints = 1; 
nOutputs = 2; % rel_mp and rel_noise
nParameterVar = 8; % 7 input parameters + 1 dummy variable

% set up min and max matrices

vary_log_range = [
    -1 1; % V
    -1 1; % k1
    -1 1; % k2
    -1 1; % kcat
    -1 1; % h_val
    -1 1; % rep_strength
    -1 1; % promoter_strength
    -1 1]; %dummy

% BASELINE VALUES
V_nom = 10^-15; % E. coli cell volume, L
k1_nom = 2*10^-4; % TetR enzyme degradation, sec^-1
k2_nom = 2*10^-4; % FFA product degradation/dilution, sec^-1
kcat_nom = 77.75; % TesA enzyme turnover rate, sec^-1
h_val_nom = 2.7; % promoter sensitivity
rep_strength_nom = 1E-3; % 1/alpha value
promoter_strength_nom = 5.86E-8; % beta value
dummy = 1;

Parameters_log_nominal = log10([V_nom, k1_nom, k2_nom, kcat_nom, h_val_nom, rep_strength_nom, promoter_strength_nom, dummy]);

pmin = vary_log_range(:,1)' + Parameters_log_nominal;
pmax = vary_log_range(:,2)' + Parameters_log_nominal;

% Parameter labels
efast_var={'V', 'k1', 'k2', 'kcat', 'h_val', 'rep_strength', 'promoter_strength', 'dummy'};

% Variable Labels
y_var_label={'rel_mp','rel_noise'};
