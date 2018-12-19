%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,GPS_DATA,...
%    HOL_RSL,HOL_AGE,HOL_LOC,priorChoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a structure filled with the hyperparameters (parameters of the
% prior distributions) following the rationale outlined in
% "Piecuch_model_description.pdf"
% INPUT: 
%   N           Number of target locations (here 211)
%   K           Number of time steps (here 118, 1900-2017)
%   M           Number of tide gauges (here 53)
%   LON         Longitudes of target locations
%   LAT         Latitudes of target locations
%   TGR_DATA    Array of tide gauge RSL time series
%   GPS_DATA    Vector of GPS vectical velocities
%   HOL_RSL     RSL values from sea level index points
%   HOL_AGE     Age values from sea level index points
%   HOL_LOC     Saltmarsh number/index for a given index point
%   priorChoice Flag to choose wide/normal/narrow priors on inverse ranges
%               (value can be 1,2,3 ... see below usage)
% OUTPUT:
%   HP          Structure of hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,GPS_DATA,...
    HOL_RSL,HOL_AGE,HOL_LOC,priorChoice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See "Piecuch_model_description.pdf" for rationale behind the priors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1:K;
m=nan(M,1);
s=nan(M,1);
r=nan(M,1);
e=nan(M,1);
y0=nan(M,1);
l=nan(M,1);
for n=1:M
    y=[]; 
    y=squeeze(TGR_DATA(n,:)); 
    i=[]; 
    i=find(~isnan(y));
        p=[]; q=[]; [p,q]=polyfit(t(i),y(i),1);
        w=(inv(q.R)*inv(q.R)')*q.normr^2/q.df; 
        m(n)=p(1);
        s(n)=w(1,1);
        [a b]=aryule(y(i)-p(1)*t(i)-p(2),1);
        d(n)=std(y(i)-p(1)*t(i)-p(2));
        r(n)=-a(2);
        e(n)=sqrt(b);
        l(n)=p(1)*mean(t)+p(2);
        y0(n)=p(2)-l(n);
    clear y i p w a b
end
u = GPS_DATA;

for nn=1:numel(unique(HOL_LOC))
 ii=[]; ii=find(HOL_LOC==nn);
 pp=[]; pp=polyfit(HOL_AGE(ii),HOL_RSL(ii),1);
 mm(nn)=pp(1);
 bb(nn)=pp(2);
 HOL_PAR(ii)=polyval(pp,HOL_AGE(ii));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variance inflation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_infl=5^2;
var_infl2=10^2;
var0_infl=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_0
% see lines 22-28 in "Piecuch_model_description.pdf"
HP.eta_tilde_y_0         	= nanmean(y0); % Mean of y_0 prior
HP.zeta_tilde_y_0_2         = var0_infl*(nanvar(y0)); % Variance of y_0 prior

% r
% see lines 48-49 in "Piecuch_model_description.pdf"
HP.u_tilde_r                = 0.1;% Lower bound of r prior
HP.v_tilde_r                = 0.9; % Upper bound of r prior

% mu
% see lines 39-43 in "Piecuch_model_description.pdf"
HP.eta_tilde_mu             = nanmean(m)+nanmean(u);% Mean of mu prior
HP.zeta_tilde_mu_2          = var_infl2*(nanvar(m)+nanvar(u)); % Variance of mu prior

% nu
% see lines 44-47 in "Piecuch_model_description.pdf"
HP.eta_tilde_nu             = nanmean(l); % Mean of nu prior
HP.zeta_tilde_nu_2          = var_infl*nanvar(l); % Variance of nu prior

% alpha
% see lines 29-32 in "Piecuch_model_description.pdf"
HP.eta_tilde_alpha          = nanmean(u); % Mean of alpha prior
HP.zeta_tilde_alpha_2       = var_infl*nanvar(u); % Variance of alpha prior

% beta
% see lines 33-38 in "Piecuch_model_description.pdf"
HP.eta_tilde_beta           = nanmean(bb); % Mean of beta prior
HP.zeta_tilde_beta_2        = var0_infl*nanvar(bb); % Variance of beta prior

% pi_2
% see lines 62-65 in "Piecuch_model_description.pdf"
HP.xi_tilde_pi_2            = 1/2; % Shape of pi_2 prior
HP.chi_tilde_pi_2           = 1/2*nanvar(m); % Inverse scale of pi_2 prior

% delta_2
% see lines 72-74 in "Piecuch_model_description.pdf"
HP.xi_tilde_delta_2         = 1/2; % Shape of delta_2 prior
HP.chi_tilde_delta_2        = 1/2*1e-4; % Guess (1 cm)^2 error variance

% sigma_2
% see lines 66-71 in "Piecuch_model_description.pdf"
HP.xi_tilde_sigma_2         = 1/2; % Shape of sigma_2 prior
HP.chi_tilde_sigma_2        = 1/2*nanmean(e.^2); % Inverse scale of sigma_2 prior

% tau_2
% see lines 85-88 in "Piecuch_model_description.pdf"
HP.xi_tilde_tau_2           = 1/2; % Shape of tau_2 prior
HP.chi_tilde_tau_2          = 1/2*nanvar(l); % Inverse scale of tau_2 prior

% gamma_2
% see lines 95-97 in "Piecuch_model_description.pdf"
HP.xi_tilde_gamma_2         = 1/2; % Shape of tau_2 prior
HP.chi_tilde_gamma_2        = 1/2*(1e-3)^2; % Guess (1 mm/yr)^2 error variance

% omega_2
% see lines 92-94 in "Piecuch_model_description.pdf"
HP.xi_tilde_omega_2         = 1/2; % Shape of tau_2 prior
HP.chi_tilde_omega_2        = 1/2*nanvar(u); % 

% epsilon_2
% see lines 89-91 in "Piecuch_model_description.pdf"
HP.xi_tilde_epsilon_2       = 1/2; % Shape of tau_2 prior
HP.chi_tilde_epsilon_2      = 1/2*(1e-3)^2; % Guess (1 mm/yr)^2 error variance

% eps_2
% see lines 81-84 in "Piecuch_model_description.pdf"
HP.xi_tilde_eps_2           = 1/2; % Shape of tau_2 prior
HP.chi_tilde_eps_2          = 1/2*1; % Guess (1 m)^2 representation error variance

% kappa_2
% see lines 75-77 in "Piecuch_model_description.pdf"
HP.xi_tilde_kappa_2         = 1/2; % Shape of tau_2 prior
HP.chi_tilde_kappa_2        = 1/10*1/2*nanvar(bb); % 

priorRangeMean              = -6.9;
% priorChoice = 1 --> (0.35)^2 -- middle of the road
% priorChoice = 2 --> (0.70)^2 -- wide priors
% priorChoice = 3 --> (0.10)^2 -- narrow priors

if priorChoice==1
    priorRangeVar = (0.35)^2;
elseif priorChoice==2
    priorRangeVar = (0.70)^2;
elseif priorChoice==3
    priorRangeVar = (0.10)^2;
else
    error('Invalid choice for priorChoice input.  Must be 1,2,3')
end

% phi
% see lines 50-53 in "Piecuch_model_description.pdf"
HP.eta_tilde_phi            = priorRangeMean; % "Mean" of phi prior
HP.zeta_tilde_phi_2         = priorRangeVar; % "Variance" of phi prior

% lambda (this one's strongly constrained; 95% within 500,2000 km)
% see lines 54-57 in "Piecuch_model_description.pdf"
HP.eta_tilde_lambda         = priorRangeMean; % "Mean" of phi prior
HP.zeta_tilde_lambda_2      = priorRangeVar; % "Variance" of phi prior

% rho (this one's strongly constrained; 95% within 500,2000 km)
% see lines 58-61 in "Piecuch_model_description.pdf"
HP.eta_tilde_rho            = priorRangeMean; % "Mean" of phi prior
HP.zeta_tilde_rho_2         = priorRangeVar; % "Variance" of phi prior

% GIA
% these priors (mean vectors and covariance matrices of GIA-driven SSH and
% VLM) were determined based on 216 GIA model predictions using different
% lithospheric thicknesses, mantle viscosities, and ice histories, as
% described in Piecuch et al. (2018).  The values loaded here correspond
% solely to the 211 target locations along the U.S. East Coast considered
% in Piecuch et al. (2018).  Readers interested in regularly gridded GIA
% model solutions should contact Jerry X. Mitrovica at Harvard University
% (jxm@eps.harvard.edu)
load GIA/GIA_priors.mat

HP.eta_tilde_ug = mu_legendreVlm;
HP.delta_tilde_ug_2 = sig_legendreVlm;
HP.inv_delta_tilde_ug_2 = inv(HP.delta_tilde_ug_2);

HP.eta_tilde_wg = mu_legendreSsh;
HP.delta_tilde_wg_2 = sig_legendreSsh;
HP.inv_delta_tilde_wg_2 = inv(HP.delta_tilde_wg_2);
