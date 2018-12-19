%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init_vals_pickup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and initialize values in case of pickup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn_nonnan = find(~isnan(R),1,'last');

mu=MU(nn_nonnan);
nu=NU(nn_nonnan);
pi_2=PI_2(nn_nonnan);
delta_2=DELTA_2(nn_nonnan);
sigma_2=SIGMA_2(nn_nonnan);
gamma_2=GAMMA_2(nn_nonnan);
tau_2=TAU_2(nn_nonnan);
phi=PHI(nn_nonnan);
lambda=LAMBDA(nn_nonnan);
alpha=ALPHA(nn_nonnan);
omega_2=OMEGA_2(nn_nonnan);
epsilon_2=EPSILON_2(nn_nonnan);
rho=RHO(nn_nonnan);
eps_2=EPS_2(nn_nonnan);
kappa_2=KAPPA_2(nn_nonnan);
beta=BETA(nn_nonnan);

a=squeeze(A(nn_nonnan,:))';
b=squeeze(B(nn_nonnan,:))';
l=squeeze(ELL(nn_nonnan,:))';
r=squeeze(R(nn_nonnan))';
y_0=squeeze(Y_0(nn_nonnan,:))';
u=squeeze(U(nn_nonnan,:))';
ug=squeeze(UG(nn_nonnan,:))';
wg=squeeze(WG(nn_nonnan,:))';
v=squeeze(V(nn_nonnan,:))';
iota=squeeze(IOTA(nn_nonnan,:))';
wye=squeeze(WYE(nn_nonnan,:))';
tee=squeeze(TEE(nn_nonnan,:))';

y=squeeze(Y(nn_nonnan,:,:))';

