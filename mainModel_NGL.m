% First order and total effect indices for a given model computed with
% Extended Fourier Amplitude Sensitivity Test (EFAST). 
% Andrea Saltelli, Stefano Tarantola and Karen Chan. % 1999.
% "A quantitative model-independent method for global sensitivity analysis
% of model output". % Technometrics 41:39-56.
clear; 
close all; 
clc; 
%% INPUT
NR = 7; %: no. of search curves - RESAMPLING
k = 8 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 600; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)

Parameter_settings_NGL;

% Computation of the frequency for the group of interest OMi and the # of
% sample points NS (here N=NS)
% Minimum NS for eFAST is 65.
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end

%% Pre-allocation of the output matrix Y
% Y will save only the points of interest specified in
% the vector time_points
Y(NS,nTimePoints,nOutputs,nParameterVar,NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for L=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        % Do the NS model evaluations.
        for run_num=1:NS
            if run_num == 1
                [i run_num L] % keeps track of [parameter run NR]
            end
            % ODE system file
            f=@modelODE_NGL;
            modelOutput=f(X(:,:,i,L),run_num);
            % ODE solver call
            % [t,y]=ode15s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,[]); 
            % It saves only the output at the time points of interest
            Y(run_num,:,:,i,L)=modelOutput;
        end %run_num=1:NS
    end % L=1:NR
end % i=1:k

save mainModel_NGL_results.mat;

% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,nTimePoints,1:nOutputs);
% Calculate Coeff. of Var. for Si and STi for rel_noise (variable 2). See
% online Supplement A.5 for details.
[CVsi CVsti]=CVmethod(Si,rangeSi,Sti,rangeSti,2);

% T-test on Si and STi for both output variables
% For rel_mp (variable 1)
s_model_1 = efast_ttest(Si,rangeSi,Sti,rangeSti,1:nTimePoints,efast_var,1,y_var_label,0.05);
% For rel_noise (variable 2)
s_model_2 = efast_ttest(Si,rangeSi,Sti,rangeSti,1:nTimePoints,efast_var,2,y_var_label,0.05);

% Calculate mean and standard deviation of Sti values across all search
% curves for both output variables. 
for i=1:k
    % Mean and STD for rel_mp (variable 1)
    meanSti_var1(i) = mean(s_model_1.rangeSti(i,:,:,1));
    stdSti_var1(i) = std(s_model_1.rangeSti(i,:,:,1));
    % Mean and STD for rel_noise (variable 2)
    meanSti_var2(i) = mean(s_model_2.rangeSti(i,:,:,2));
    stdSti_var2(i) = std(s_model_2.rangeSti(i,:,:,2));
end

%% FUNCTIONS
% Function defined in SETFREQ.m
% Algorithm for selection of a frequency set for the complementary group.
% Done recursively as described in: Appendix of "Sensitivity Analysis"
% [Saltelli et al., 2000]
% OMci = SETFREQ(k-1,OMi/2/MI)
function OMci = SETFREQ(Kci,OMciMAX,f)
if Kci==1
    OMci = 1;
elseif OMciMAX==1
    OMci = ones(1,Kci);
else
    if(OMciMAX < Kci)
         INFD = OMciMAX;
    else
        INFD = Kci;
    end
    ISTEP = round((OMciMAX-1)/(INFD-1));
    if(OMciMAX == 1)
        ISTEP = 0;
    end
    OTMP = 1:ISTEP:INFD*ISTEP;
    fl_INFD = floor(INFD);
    for i=1:Kci
        j = mod(i-1,fl_INFD)+1;
        OMci(i) = OTMP(j);
    end
end
OMci(f)=[];
end

% Function defined in CVmethod.m
function [CVsi CVsti]=CVmethod(Si,rangeSi,Sti,rangeSti,out)
meanSi=[];
meanSti=[];
[k s NR u]=size(rangeSi);
u;
if u==1
    out=1;
    for j=1:k
        for t=1:s
            meanSi(j,t,out)=(mean(rangeSi(j,t,:,out)));
            meanSti(j,t,out)=(mean(rangeSti(j,t,:,out)));
            stdSi(j,t,out)=(std(rangeSi(j,t,:,out)));
            stdSti(j,t,out)=(std(rangeSti(j,t,:,out)));
        end
    end
    a=Si(:,:,out)./meanSi(:,:,out);
    b=Sti(:,:,out)./meanSti(:,:,out);
else
    for j=1:k
        for t=1:s
            meanSi(j,t,out)=squeeze(mean(rangeSi(j,t,:,out)));
            meanSti(j,t,out)=squeeze(mean(rangeSti(j,t,:,out)));
            stdSi(j,t,out)=squeeze(std(rangeSi(j,t,:,out)));
            stdSti(j,t,out)=squeeze(std(rangeSti(j,t,:,out)));
        end
    end
a=stdSi(:,:,out)./meanSi(:,:,out);
b=stdSti(:,:,out)./meanSti(:,:,out);
    %a=Si(:,:,out)./meanSi(:,:,out)
    %b=Sti(:,:,out)./meanSti(:,:,out)
end
%for i=1:k
%    for t=1:s
%        CVsi(i,t)=std(a(i,t))/mean(a(i,t));
%        CVsti(i,t)=std(b(i,t))/mean(b(i,t));
%    end
%end

CVsi=100*a';
CVsti=100*b';
end

% Function defined in parameterdist.m
% Parameter distributions for eFAST sampling scheme %%
% X: efast search curves normalized between o and 1, i.e. 
% [X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi;
% [pmin pmax]: min and max values of the range of variation
% [pmean pstd]: mean and standard deviation for distributions other than
% uniform
% nsample: number of samples (usually 65 in efast)
% type: type of distributions (pdf). Uniform ['unif'], nomral ['norm'],
% lognormal ['lognorm']
% The uniform pdf implements a log scale for pmax/pim>1e4
function Xdist = parameterdist(X,pmax,pmin,pmean,pstd,nsample,type)
switch lower(type)
    case {'unif'}
        for k=1:length(X(1,:)) %loop through parameters
            nvar=length(pmin);
            nsample;
            ran=rand(nsample,1);
            s=zeros(nsample,1);
            idx=randperm(nsample);
            %pause;
            P =(idx'-ran)/nsample;
            %pause;
            Xdist(:,k)=(X(:,k).*(pmax(k)-pmin(k)))+pmin(k);
        end
    case {'norm'}
        for k=1:length(X(1,:)) %loop through parameters
            Xdist(:,k) = norminv(X(:,k),pmean(k),pstd(k));
        end
    case {'lognorm'}
        for k=1:length(X(1,:)) %loop through parameters
            Xdist(:,k) = norminv(X(:,k),log(pmean(k)),pstd(k));
        end
    otherwise
        disp('Unknown pdf')
end
end