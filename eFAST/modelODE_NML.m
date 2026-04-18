% This ODE represents the negative feedback loop with indirect repression 
% of additional TesA expression via FFA binding of FadR to repress
% further activation of the hybrid promoter.
function [modelOutput] = modelODE_NML(P,run_num)

% PARAMETER INITIALIZATION
V = 10.^P(run_num,1); % E. coli cell volume, L
k1 = 10.^P(run_num,2); % TesA enzyme degradation/dilution, sec^-1
k2 = 10.^P(run_num,3); % FFA product degradation/dilution, sec^-1
kcat = 10.^P(run_num,4); % TesA enzyme turnover rate, sec^-1
h_val = 10.^P(run_num,5); % promoter sensitivity
rep_strength = 10.^P(run_num,6); % 1/alpha values
promoter_strength = 10.^P(run_num,7); % beta values 

Na = 6.022*10^23; % Avogadro's number, molecules mol^-1
k = k1*k2/kcat;
alpha = 1/rep_strength;
beta = promoter_strength;
h = h_val;

% syms mu mp

% Calculate mean product number
SSE_prod = @(mu)(k*mu+k*(mu^(h+1))/(alpha^h)-beta)^2;
mean_prod = fminsearch(SSE_prod, 0.1);
rel_mp = mean_prod*k/beta;

% Calculate relative product noise
mu = mean_prod;
rel_noise = calculateRelNoise(h, beta, alpha, mu, k1, k2, kcat);

% Final output is a single value for rel_mp and rel_noise
modelOutput = [rel_mp rel_noise];

% FUNCTIONS
function y = calculateRelNoise(h, beta, alpha, mu, k1, k2, kcat)
    k = k1*k2/kcat;
    y = (beta/(k*mu))-((k1+kcat)/(k1+k2+kcat))*(beta*h*(mu^h)/(beta*(alpha^h)+h*k*(mu^(h+1))));
end

end
