%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize vectors and matrices of posterior Bayes model solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = nan(NN_thin,M);         % Spatial field of error trend
B = nan(NN_thin,N);         % Spatial field of process linear trend
ELL = nan(NN_thin,M);       % Spatial field of observational biases
R = nan(NN_thin,1);         % AR(1) coefficient of the process
Y = nan(NN_thin,K,N);       % Process values
Y_0 = nan(NN_thin,N);       % Process initial conditions
MU = nan(NN_thin,1);        % Mean value of process linear trend
NU = nan(NN_thin,1);        % Mean value of observational biases
PHI = nan(NN_thin,1);       % Inverse range of process innovations
LAMBDA = nan(NN_thin,1);    % Inverse range of process trend field
PI_2 = nan(NN_thin,1);      % Spatial variance of process linear trend
SIGMA_2 = nan(NN_thin,1);   % Sill of the process innovations
DELTA_2 = nan(NN_thin,1);   % Instrumental error variance 
TAU_2 = nan(NN_thin,1);     % Spatial variance in observational biases
GAMMA_2 = nan(NN_thin,1);   % Spatial variance in error trends
KAPPA_2 = nan(NN_thin,1);   % Spatial variance in index-point intercepts
EPS_2 = nan(NN_thin,1);     % Spatial variance in index-point residuals

BETA = nan(NN_thin,1);      % Spatial mean in index-point intercepts
IOTA = nan(NN_thin,N_S);    % Vector of index-point intercepts

WYE = nan(NN_thin,N_D);     % Vector of index-point process values
TEE = nan(NN_thin,N_D);     % Vector of index-point times

U = nan(NN_thin,N);         % Vector of regional VLM rate
UG= nan(NN_thin,N);         % Vector of GIA-driven VLM rate
WG= nan(NN_thin,N);         % Vector of GIA-driven SSH rate
V = nan(NN_thin,N);         % Total VLM vector
ALPHA = nan(NN_thin,1);     % Mean non-GIA VLM rate
OMEGA_2 = nan(NN_thin,1);   % Partial sill of regional VLM rate
RHO = nan(NN_thin,1);       % non-GIA VLM rate inverse range
EPSILON_2 = nan(NN_thin,1); % Variance in tide gauge data biases
