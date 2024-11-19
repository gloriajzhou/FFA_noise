% ngl_solv contains the code for modeling the NGL system using parameters
% specific to FFA synthesis in Escherichia coli in minimal glucose. The NGL
% system is a direct form of feedback repression from polycistronic TetR
% that reduces further TesA expression. 

clear;
close all;
clc;

% Mean enzyme number (fixed h value)
Na = 6.022*10^23; % Avogadro's number, molecules mol^-1
V = 10^-15; % E. coli cell volume, L
k1 = 4*10^-4; % TetR enzyme degradation, sec^-1
k2 = 2*10^-4; % FFA product degradation/dilution, sec^-1
kcat = 77.75; % TesA enzyme turnover rate, sec^-1
k = k1*k2/kcat;

h_val = 2.7; % promoter sensitivity
n_val = 1; % TesA/TetR scaling factor 
rep_strength = logspace(-4, -2, 10); % 1/alpha values
promoter_strength = logspace(-3, 2, 100); % beta values
mean_enz = []; % initialize empty vector
rel_mean = [];
syms mu

for i = 1:length(rep_strength)
    alpha = 1/rep_strength(i);
    for j = 1:length(promoter_strength)
        beta = promoter_strength(j); % molecules sec^-1
        h = h_val;
        SSE = @(mu)(k1*mu+(k1*(mu^(h+1))/(alpha^h))-beta)^2;
        mean_enz(i, j) = fminsearch(SSE, 0.1);
        % rel_mean(i,j) = mean_enz(i,j)*(k1/beta);
    end
end

% figure
% plot(promoter_strength, mean_enz)
% xlabel('Promoter strength');
% ylabel('Mean enzyme number');

% Mean product number (fixed h value)
mean_prod = []; % initialize empty vector
rel_mp = [];
syms mp

for i = 1:length(rep_strength)
    alpha = 1/rep_strength(i);
    for j = 1:length(promoter_strength)
        beta = promoter_strength(j);
        h = h_val;
        n = n_val;
        SSE = @(mp)(k*mp/n+(k1/(alpha^h))*((k2*mp/(kcat*n))^(h+1))-beta)^2;
        mean_prod(i, j) = fminsearch(SSE, 0.1);
        rel_mp(i,j) = mean_prod(i,j)*k/(beta*n);
    end
end

% for i = 1:length(promoter_strength)
%     alpha = 10^1;
%     beta = promoter_strength(i);
%     h = h_val;
%     n = n_val;
%     SSE = @(mp)(k*mp/n+(k5/(alpha^h))*((k6*mp/(kcat*n))^(h+1))-beta)^2;
%     mean_prod(i) = fminsearch(SSE, 0.1);
%     rel_mp(i) = mean_prod(i)*k/(beta*n);
% end
% 
% figure
% loglog(promoter_strength, rel_mp)
% xlabel('Promoter strength');
% ylabel('Relative mean product number');

% Relative product noise
rel_noise = []; % initialize empty vector

for i = 1:length(rep_strength)
    alpha = 1/rep_strength(i);
    for j = 1:length(promoter_strength)
        mu = mean_enz(i, j);
        beta = promoter_strength(j);
        h = h_val;
        n = n_val;
        value = calculateRelNoise(h, beta, alpha, mu, n);
        rel_noise(i, j) = value;
    end
end

% Create figures for rel mean and rel noise
figure 
winter_colormap = winter(length(rep_strength));
% yyaxis left

for i = 1:length(rep_strength)
    color = winter_colormap(i, :);
    loglog(promoter_strength, rel_noise(i, :), 'Color', color, 'LineWidth', 2, 'LineStyle','-','Marker', 'none');
    hold on;
end

colormap(winter);
c = colorbar;
ylabel(c, 'Repression strength [1/molecules]');
clim([min(rep_strength) max(rep_strength)]);
xlabel('Promoter strength [molecules/sec]');
ylabel('Relative product noise')

% Figure for relative mean product number
figure 
winter_colormap = winter(length(rep_strength));

% yyaxis right
% hold on 

for i = 1:length(rep_strength)
    color = winter_colormap(i, :);
    loglog(promoter_strength, rel_mp(i, :), 'Color', color, 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none');
    hold on;
end

colormap(winter);
c = colorbar;
ylabel(c, 'Repression strength [1/molecules]');
clim([min(rep_strength) max(rep_strength)]);

% title('Fixed h value');
xlabel('Promoter strength [molecules/sec]');
ylabel('Relative mean product number');

hold off

% Define function
function y = calculateRelNoise(h, beta, alpha, mu, n)
    k1 = 4*10^-4; % enzyme degradation, sec^-1
    k2 = 2*10^-4; % product consumption, sec^-1
    kcat = 77.75; % enzyme turnover rate, sec^-1    
    fderiv = h*(k1^2)*(mu^(h+1))/(beta*(alpha^h));
    y = (beta/(k1*mu))*(1+(kcat*n*k1)/((k1+k2+fderiv)*(k1+fderiv)))/(1+(kcat*n/(k1+k2)));
end

