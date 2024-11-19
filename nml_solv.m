% nml_solv contains the code for modeling the NML system using parameters
% specific to FFA synthesis in Escherichia coli in minimal glucose. The NML
% system is a product-sensitive approach to repressing further TesA
% expression through inhibition of further FadR activity by FFA.

% Mean product number (fixed h value)
Na = 6.022*10^23; % Avogadro's number, molecules mol^-1
V = 10^-15; % E. coli cell volume, L
k1 = 2*10^-4; % enzyme degradation, sec^-1
k2 = 2*10^-4; % product degradation/dilution, sec^-1
kcat = 77.75; % TesA enzyme turnover rate, sec^-1
k = k1*k2/kcat;

h_val = 2.5; % promoter sensitivity
rep_strength = logspace(-4, -2, 10); % 1/alpha values
promoter_strength = logspace(-9, -4, 100); % beta values

mean_prod = []; % initialize empty vector
rel_mean = [];
syms mu

for i = 1:length(rep_strength)
    alpha = 1/rep_strength(i);
    for j = 1:length(promoter_strength)
        beta = promoter_strength(j);
        h = h_val;
        SSE = @(mu)(k*mu+k*(mu^(h+1))/(alpha^h)-beta)^2;
        mean_prod(i, j) = fminsearch(SSE, 0.1);
        rel_mean(i,j) = mean_prod(i,j)*(k/beta);
    end
end

% Relative noise (fixed h value)
rel_noise = []; % initialize empty vector

for i = 1:length(rep_strength)
    alpha = 1/rep_strength(i);
    for j = 1:length(promoter_strength)
        h = h_val;
        beta = promoter_strength(j);
        mu = mean_prod(i, j);
        value = calculateRelNoise(h, beta, alpha, mu);
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
    loglog(promoter_strength, rel_mean(i, :), 'Color', color, 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none');
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

% Define functions
function y = calculateRelNoise(h, beta, alpha, mu)
    k1 = 2*10^-4; % enzyme degradation, sec^-1
    k2 = 2*10^-4; % product consumption, sec^-1
    kcat = 77.75; % enzyme turnover rate, sec^-1
    k = k1*k2/kcat;
    y = (beta/(k*mu))-((k1+kcat)/(k1+k2+kcat))*(beta*h*(mu^h)/(beta*(alpha^h)+h*k*(mu^(h+1))));
end
