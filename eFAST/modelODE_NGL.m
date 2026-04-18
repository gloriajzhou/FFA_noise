% This ODE represents the negative feedback loop with direct repression of
% additional TesA expression via TetR binding to the hybrid promoter.
function [modelOutput] = modelODE_NGL(P,run_num)

% PARAMETER INITIALIZATION
V = 10.^P(run_num,1); % E. coli cell volume, L
k1 = 10.^P(run_num,2); % TetR enzyme degradation, sec^-1
k2 = 10.^P(run_num,3); % FFA product degradation/dilution, sec^-1
kcat = 10.^P(run_num,4); % TesA enzyme turnover rate, sec^-1
h_val = 10.^P(run_num,5); % promoter sensitivity
n_val = 10.^P(run_num,6); % TesA/TetR scaling factor 
rep_strength = 10.^P(run_num,7); % 1/alpha values
promoter_strength = 10.^P(run_num,8); % beta values 

Na = 6.022*10^23; % Avogadro's number, molecules mol^-1
k = k1*k2/kcat;
alpha = 1/rep_strength;
beta = promoter_strength;
h = h_val;
n = n_val;
% syms mu mp

% Calculate mean enzyme number
SSE_enz = @(mu)(k1*mu+(k1*(mu^(h+1))/(alpha^h))-beta)^2;
mean_enz = fminsearch(SSE_enz, 0.1);

% Calculate mean product number
SSE_prod = @(mp)(k*mp/n+(k1/(alpha^h))*((k2*mp/(kcat*n))^(h+1))-beta)^2;
mean_prod = fminsearch(SSE_prod, 0.1);
rel_mp = mean_prod*k/(beta*n);

% Calculate relative product noise
mu = mean_enz;
rel_noise = calculateRelNoise(h, beta, alpha, mu, n, k1, k2, kcat);

% Final output is a single value for rel_mp and rel_noise
modelOutput = [rel_mp rel_noise];

% FUNCTIONS
function y = calculateRelNoise(h, beta, alpha, mu, n, k1, k2, kcat)
    fderiv = h*(k1^2)*(mu^(h+1))/(beta*(alpha^h));
    y = (beta/(k1*mu))*(1+(kcat*n*k1)/((k1+k2+fderiv)*(k1+fderiv)))/(1+(kcat*n/(k1+k2)));
end

end
