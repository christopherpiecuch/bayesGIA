%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete_burn_in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deletes specified warm-up interations from solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BETA(1:NN_burn_thin)=[];
EPS_2(1:NN_burn_thin)=[];
KAPPA_2(1:NN_burn_thin)=[];
IOTA(1:NN_burn_thin,:)=[];
TEE(1:NN_burn_thin,:)=[];
WYE(1:NN_burn_thin,:)=[];


MU(1:NN_burn_thin)=[];
NU(1:NN_burn_thin)=[];
PI_2(1:NN_burn_thin)=[];
DELTA_2(1:NN_burn_thin)=[];
SIGMA_2(1:NN_burn_thin)=[];
TAU_2(1:NN_burn_thin)=[];
GAMMA_2(1:NN_burn_thin)=[];
LAMBDA(1:NN_burn_thin)=[];
PHI(1:NN_burn_thin)=[];
A(1:NN_burn_thin,:)=[];
B(1:NN_burn_thin,:)=[];
ELL(1:NN_burn_thin,:)=[];
R(1:NN_burn_thin)=[];
Y_0(1:NN_burn_thin,:)=[];
Y(1:NN_burn_thin,:,:)=[];

U(1:NN_burn_thin,:)=[];
UG(1:NN_burn_thin,:)=[];
WG(1:NN_burn_thin,:)=[];
V(1:NN_burn_thin,:)=[];
ALPHA(1:NN_burn_thin)=[];
RHO(1:NN_burn_thin)=[];
OMEGA_2(1:NN_burn_thin)=[];
EPSILON_2(1:NN_burn_thin)=[];
