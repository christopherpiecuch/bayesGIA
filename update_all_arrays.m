%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_all_arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store latest value of parameter and process values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BETA(nn_thin)=beta;
KAPPA_2(nn_thin)=kappa_2;
EPS_2(nn_thin)=eps_2;
MU(nn_thin)=mu;
NU(nn_thin)=nu;
PI_2(nn_thin)=pi_2;
DELTA_2(nn_thin)=delta_2;
SIGMA_2(nn_thin)=sigma_2;
GAMMA_2(nn_thin)=gamma_2;
TAU_2(nn_thin)=tau_2;
PHI(nn_thin)=phi;
LAMBDA(nn_thin)=lambda;

IOTA(nn_thin,:)=iota;
WYE(nn_thin,:)=wye;
TEE(nn_thin,:)=tee;


A(nn_thin,:)=a;
B(nn_thin,:)=b;
ELL(nn_thin,:)=l;
R(nn_thin)=r;
Y_0(nn_thin,:)=y_0;
Y(nn_thin,:,:)=y';

U(nn_thin,:)=u;
UG(nn_thin,:)=ug;
WG(nn_thin,:)=wg;
V(nn_thin,:)=v;
ALPHA(nn_thin)=alpha;
OMEGA_2(nn_thin)=omega_2;
EPSILON_2(nn_thin)=epsilon_2;
RHO(nn_thin)=rho;