%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set_initial_values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define initial process and parameter values as described in methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v

% mean parameters
mu=normrnd(HP.eta_tilde_mu,sqrt(HP.zeta_tilde_mu_2));
nu=normrnd(HP.eta_tilde_nu,sqrt(HP.zeta_tilde_nu_2));
alpha=[];
alpha=normrnd(HP.eta_tilde_alpha,sqrt(HP.zeta_tilde_alpha_2));
beta=normrnd(HP.eta_tilde_beta,sqrt(HP.zeta_tilde_beta_2));

% variance parameters
pi_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_pi_2,HP.xi_tilde_pi_2], [1,1])]); % use min to prevent needlessly large values
delta_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_delta_2,HP.xi_tilde_delta_2], [1,1])]); % use min to prevent needlessly large values
sigma_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_sigma_2,HP.xi_tilde_sigma_2], [1,1])]); % use min to prevent needlessly large values
tau_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_tau_2,HP.xi_tilde_tau_2], [1,1])]); % use min to prevent needlessly large values
gamma_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_gamma_2,HP.xi_tilde_gamma_2], [1,1])]); % use min to prevent needlessly large values
omega_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_omega_2,HP.xi_tilde_omega_2], [1,1])]); % use min to prevent needlessly large values
epsilon_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_epsilon_2,HP.xi_tilde_epsilon_2], [1,1])]); % use min to prevent needlessly large values
eps_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_eps_2,HP.xi_tilde_eps_2], [1,1])]); % use min to prevent needlessly large values
kappa_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_kappa_2,HP.xi_tilde_kappa_2], [1,1])]); % use min to prevent needlessly large values

% inverse length scale parameters
phi=exp(normrnd(HP.eta_tilde_phi,sqrt(HP.zeta_tilde_phi_2)));
lambda=exp(normrnd(HP.eta_tilde_lambda,sqrt(HP.zeta_tilde_lambda_2)));
rho=exp(normrnd(HP.eta_tilde_rho,sqrt(HP.zeta_tilde_rho_2)));

% spatial fields
a=(mvnrnd(zeros(M,1),gamma_2*eye(M)))';        
iota=(mvnrnd(beta*ones(N_S,1),kappa_2*eye(N_S)))';        
b=(mvnrnd(mu*ones(N,1),pi_2*exp(-lambda*D)))';        
l=(mvnrnd(nu*ones(M,1),tau_2*eye(M)))';        
u=(mvnrnd(alpha*ones(N,1),omega_2*exp(-rho*D)))';        
v=(mvnrnd(u,epsilon_2*eye(N)))';        
ug=(mvnrnd(HP.eta_tilde_ug,HP.delta_tilde_ug_2))';        
wg=(mvnrnd(HP.eta_tilde_wg,HP.delta_tilde_wg_2))';        

tee=(mvnrnd(S,XiMat))';
diagtee=diag(tee);
wye=(mvnrnd(H,GammaMat))';        

% AR(1) parameter
r=HP.u_tilde_r+(HP.v_tilde_r-HP.u_tilde_r )*rand(1);

% process
y_0=zeros(N,1);
y=zeros(N,K);
