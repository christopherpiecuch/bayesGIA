%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_main_code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main driver MATLAB code for producing results in the main text of
%   Piecuch, Huybers, Hay, Kemp, Little, Mitrovica, Ponte, and Tingley,
%   "Origin of spatial variation in United States East Coast sea level
%   trends during 1900-2017", 2018, Nature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_main_code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Say hello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(0.1)
disp('Hello.  Things have started.')
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some preliminary input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% name of experiment and file to be saved
nameOfExperiment='test';                % experiment name
save_name=[date,'_',nameOfExperiment];  % file name

%%%
% number of draws to perform
NN_burn=200;             % 200,000 warm-up draws
NN_post=200;             % 200,000 post-warm-up draws
thin_period=1;            % thin chains keeping 1 of 200
startFromPickup=0;          % start from existing pickup or not
                            % default=0
%%%
% geographic bounds
inp.la1=24;    				% Southern latitudinal bounds of study region
inp.la2=46;    				% Northern latitudinal bounds "
inp.lo1=-81.8067;    		% Western longitudinal bounds "
inp.lo2=-66.5;    			% Eastern longitudinal bounds "

%%%
% minimum length of tide gauge record to be involved in the estimation
inp.minnum=25;  			% units are years

%%%
% PSMSL coastline codes to be considered (see http://www.psmsl.org/data/obtaining/)
inp.coastcode=[940 960]; 		% PSMSL ID for US east and gulf coast
                     
%%%
% PSMSL ID of additional tide gauges (not falling along the coastlines 
% specified, outside of the specified geographic bounds above, or not 
% meeting the minimum number of years requirement) to include in the 
% estimation.  See online-only methods for the rationale.  See the PSMSL
% website (see http://www.psmsl.org/data/obtaining/) for the station
% records corresponding to these ID numbers
inp.exceptions=[270 1182 2123 1696 96 1158 195 520 1106 1107]; 

%%%
% start year of inference
inp.t0 = 1900;

%%%
% prior on inverse range parameters
inp.priorChoice=1;			% prior range of ~500,2000 km (default)
%%inp.priorChoice=2;		% prior range of ~250,4000 km (wide priors)
%%inp.priorChoice=3;		% prior range of ~800,1200 km (narrow priors)

%%%
% perform joint or decoupled inference; see supplementary materials
inp.coupled=1;				% coupled inference (default)
%%inp.coupled=0; 			% decoupled inference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define iteration parameters based on input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN_burn_thin=NN_burn/thin_period; % Total number of burn-in to keep
NN_post_thin=NN_post/thin_period; % Total number of post-burn-in to keep
NN=NN_burn+NN_post;               % Total number of draws to take 
NN_thin=NN_burn_thin+NN_post_thin;% Total number of draws to keep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This line of code is for reproducing the
% coupled versus decoupled comparison in the
% supplementary information.  Default for 
% reproducing results in the main text is to
% us the ``coupled'' algorithm (calg=1).  
% Choose calg=0 to decouple the algorithm.
% See supplementary material of Piecuch 
% et al. (2018) for more information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if inp.coupled==1
    calg=1;
elseif inp.coupled==0
    calg=0;
else
    error('inp.coupled equals 0 or 1') 
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare and load PSMSL annual tide gauge RSL records, GPS VLM 
% rate data, and proxy RSL reconstructions from radiocarbon-dated
% sediment from saltmarshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% load tide gauge data
[TGR_DATA,TGR_LON,TGR_LAT,TGR_NAME] = prepare_tgr_data(inp.la1,...
    inp.la2,inp.lo1,inp.lo2,inp.minnum,inp.coastcode,inp.exceptions,inp.t0);

%%%
% load GPS data    
[GPS_DATA,GPS_ERROR,GPS_LON,GPS_LAT,GPS_NAME] = prepare_gps_data(...
    inp.la1,inp.la2,inp.lo1,inp.lo2);

%%%
% include 0.5-degree regular grid cells within study domain to estimate
% the solution at
[G_LON,G_LAT] = grab_coastal_gia_cells(inp.la1,inp.la2,inp.lo1,inp.lo2,...
    TGR_LON,TGR_LAT,GPS_LON,GPS_LAT,0.5);

%%%
% load sea level index poitns from saltmarsh sediment
load('SLIPs/late_holocene_index_points_for_bayes_2000_150_yrBP.mat')
ii=[]; ii=find(HOL_LAT>=inp.la1&HOL_LAT<=inp.la2&HOL_LON>=inp.lo1&HOL_LON<=inp.lo2);
 
% grab lats, lons, names
HOL_LON=HOL_LON(ii);
HOL_LAT=HOL_LAT(ii);
HOL_NAME=HOL_NAME(ii);
% grab data
ii=find(ismember(HOL_LOC,ii));
HOL_AGE=HOL_AGE(ii);
HOL_AGE_ERR=HOL_AGE_ERR(ii);
HOL_LOC=HOL_LOC(ii);
HOL_RSL=HOL_RSL(ii);
HOL_RSL_ERR=HOL_RSL_ERR(ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define space and time parameters related to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up dimensions of process
[N_tgr,K_tgr]=size(TGR_DATA);
[N_gps,K_gps]=size(GPS_DATA');
[N_hol,K_hol]=size(HOL_LAT);
[N_gia,K_gia]=size(G_LON);

% bring together all target locations
N = N_tgr+N_gps+N_gia+N_hol;
K = K_tgr;
LON = [TGR_LON GPS_LON G_LON' HOL_LON'];
LAT = [TGR_LAT GPS_LAT G_LAT' HOL_LAT'];

% define some useful space and time values
N_D = numel(HOL_RSL);               % number of saltmarshes
N_S = numel(HOL_LAT);               % number of RSL index points
D=EarthDistances([LON' LAT']);      % distances between locations
T=1:K; T=T-mean(T);                 % time steps (centered on zero)
T0=T(1)-1;                          % initial time step (relative to centered series)
M=sum(sum(~isnan(TGR_DATA'))~=0);   % number of locations with at least one datum
L = sum(~isnan(GPS_DATA));          % number of crustal motion data sites
% define names of all target locations
for n=1:N
    if n<=N_tgr
        NAME(n).name = TGR_NAME(n).name;
    elseif n>N_tgr&&n<=(N_tgr+N_gps)
        NAME(n).name = GPS_NAME(n-N_tgr,:);       
    elseif n>(N_tgr+N_gps)&&n<=(N_tgr+N_gps+N_gia)
        NAME(n).name = [num2str(abs(G_LON(n-N_tgr-N_gps))),'W_',...
            num2str(abs(G_LAT(n-N_tgr-N_gps))),'N'];
    else % n>(N_tgr+N_gps+N_gia)&&n<=(N_tgr+N_gps+N_gia+N_hol)
        NAME(n).name = HOL_NAME(n-N_tgr-N_gps-N_gia).name;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster locations north and south of Cape Hatteras (see Equation 
% 2 in online-only methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = determine_clusters(LAT);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the seeds of the random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(sum(clock))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate space for the sample arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=[];
initialize_output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HP = set_hyperparameters(N,K,M,LON,LAT,TGR_DATA,GPS_DATA,HOL_RSL,...
    HOL_AGE,HOL_LOC,inp.priorChoice);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename holocene index points and 
% define error covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = HOL_RSL;
S = HOL_AGE;
GammaMat=diag(HOL_RSL_ERR.^2);
XiMat=diag(HOL_AGE_ERR.^2);
invGammaMat=inv(GammaMat);
invXiMat=inv(XiMat);

%%%%%%%%%%%%%%%%%%%%
% Set initial values
%%%%%%%%%%%%%%%%%%%%
set_initial_values

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up selection matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_master=double(~isnan(TGR_DATA));
M_k=sum(H_master);
E = zeros(L,N);
E(1:L,(N_tgr+1):(N_tgr+L))=eye(L);
x = GPS_DATA';
DeltaMat = diag(GPS_ERROR.^2);
invDeltaMat = inv(DeltaMat);
for k=1:K
    gauges_with_data(k).indices=find(H_master(:,k)~=0);
    selMat(k).H=zeros(M_k(k),N);
    selMat(k).F=zeros(M_k(k),M);
    for m_k=1:M_k(k)
       selMat(k).H(m_k,gauges_with_data(k).indices(m_k))=1;
       selMat(k).F(m_k,gauges_with_data(k).indices(m_k))=1;
    end
    Z(k).z=squeeze(TGR_DATA(gauges_with_data(k).indices,k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up identity matrices and vectors of zeros or ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_N=eye(N);
I_M=eye(M);
ONE_N=ones(N,1);
ONE_M=ones(M,1);
ZERO_N=zeros(N,1);
ZERO_M=zeros(M,1);

I_N_S=eye(N_S);
I_N_D=eye(N_D);
ONE_N_S=ones(N_S,1);
ONE_N_D=ones(N_D,1);
ZERO_N_S=zeros(N_S,1);
ZERO_N_D=zeros(N_D,1);

for k=1:K
   I_MK(k).I=eye(M_k(k));
   ONE_MK(k).ONE=ones(M_k(k),1);
   ZERO_MK(k).ZERO=zeros(M_k(k),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up more selection matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selG=zeros(N_D,N);
selD=zeros(N_D,N_S);
for i=1:N_D
    for j=1:N
        if (HOL_LOC(i)+N_tgr+N_gps+N_gia)==j
            selG(i,j)=1;
        end
    end
end

for i=1:N_D
    for j=1:N_S
        if HOL_LOC(i)==j
            selD(i,j)=1;
        end
    end
end


nn_start=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Bayes code takes a while to run.  So, the code periodically saves out
% the solution "so far".  This code immediately below is to, in the case
% that the model crashes, start from a pickup file, so as to prevent having
% to start from scratch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if startFromPickup
	load(['bayes_model_solutions/experiment_',save_name,'.mat'],...
        'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
        'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'U','UG','WG','V','ALPHA',...
        'RHO','OMEGA_2','EPSILON_2','HP','*DATA','*LON','*LAT',...
        '*NAME','*ERROR','N','K','M','L','D','nn','HOL_*',...
        'KAPPA_2','EPS_2','BETA','IOTA','WYE','TEE','inp')
    
    % set start time
    nn_start=nn+1; clear nn
    
    % reset initial values
    init_vals_pickup
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through the Gibbs sampler with Metropolis step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

for nn=nn_start:NN, 
    if mod(nn,100)==0, 
        toc, 
        disp([num2str(nn),' of ',num2str(NN),' iterations done.']), 
        tic, 
    end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);
    
    
    if (mod(nn,(NN_post/10))==0)
        % periodically output, just in case of stall or crash
        disp('hello') 
        save(['bayes_model_solutions/experiment_',save_name,'.mat'],...
            'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
            'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'U','UG','WG','V','ALPHA',...
            'RHO','OMEGA_2','EPSILON_2','HP','*DATA','*LON','*LAT',...
            '*NAME','*ERROR','N','K','M','L','D','nn','HOL_*',...
            'KAPPA_2','EPS_2','BETA','IOTA','WYE','TEE','inp')
    end

     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define matrices to save time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SigmaMat=sigma_2*(C.*exp(-phi*D)); invSigmaMat=inv(SigmaMat);
    SMat=C.*exp(-phi*D); invSMat=inv(SMat);
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat); 
    LMat=exp(-lambda*D); invLMat=inv(LMat);
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat);
    TMat=exp(-rho*D); invTMat=inv(TMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_K|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 9-12 in Piecuch_model_description.pdf
    V_Y_K=[]; PSI_Y_K=[];
    V_Y_K=delta_2^(-1)*(selMat(K).H'*(Z(K).z-selMat(K).F*(l+a*T(K))))+...
    	invSigmaMat*(r*y(:,K-1)+(T(K)-r*T(K-1))*b);
    %PSI_Y_K=(1/delta_2*selMat(K).H'*selMat(K).H+invSigmaMat)^(-1);
    PSI_Y_K=inv(1/delta_2*selMat(K).H'*selMat(K).H+invSigmaMat);
    y(:,K)=mvnrnd(PSI_Y_K*V_Y_K,PSI_Y_K)';
    clear V_Y_K PSI_Y_K   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_k|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 5-8 in Piecuch_model_description.pdf
    for kk=(K-1):-1:1
    	V_Y_k=[]; PSI_Y_k=[];
      	if kk==1
        	V_Y_k=1/delta_2*(selMat(1).H'*(Z(1).z-selMat(1).F*(l+a*T(1))))+...
                invSigmaMat*(r*(y_0+y(:,2))+(1+r^2)*T(1)*b-r*(T0+T(2))*b);
        else
         	V_Y_k=1/delta_2*(selMat(kk).H'*(Z(kk).z-selMat(kk).F*(l+a*T(kk))))+...
            	invSigmaMat*(r*(y(:,kk-1)+y(:,kk+1))+(1+r^2)*T(kk)*b-r*(T(kk-1)+T(kk+1))*b);
        end
       	PSI_Y_k=inv(1/delta_2*selMat(kk).H'*selMat(kk).H+(1+r^2)*invSigmaMat);
       	y(:,kk)=mvnrnd(PSI_Y_k*V_Y_k,PSI_Y_k)';
      	clear V_Y_k PSI_Y_k 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 1-4 in Piecuch_model_description.pdf
    V_Y_0=[]; PSI_Y_0=[];
    V_Y_0=(HP.eta_tilde_y_0/HP.zeta_tilde_y_0_2)*ONE_N+invSigmaMat*(r*y(:,1)-r*(T(1)-r*T0)*b);
    PSI_Y_0=inv(1/HP.zeta_tilde_y_0_2*I_N+r^2*invSigmaMat);
    y_0=mvnrnd(PSI_Y_0*V_Y_0,PSI_Y_0)';
    clear V_Y_0 PSI_Y_0
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(a|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % see equations 17-20 in Piecuch_model_description.pdf
  	V_A=[]; PSI_A=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M);
    for kk=1:K
    	SUM_K1=SUM_K1+T(kk)*selMat(kk).F'*(Z(kk).z-...
         	selMat(kk).H*y(:,kk)-selMat(kk).F*l);
        SUM_K2=SUM_K2+T(kk)^2*(selMat(kk).F'*selMat(kk).F);
   	end
    V_A=(1/delta_2)*SUM_K1;
    PSI_A=inv(1/gamma_2*eye(M)+1/delta_2*SUM_K2);
    a=mvnrnd(PSI_A*V_A,PSI_A)';
    clear V_A PSI_A SUM_K*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(delta_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 90-92 in Piecuch_model_description.pdf
    SUM_K=0;
    for kk=1:K
      	xxx=[]; xxx=(Z(kk).z-selMat(kk).H*y(:,kk)-selMat(kk).F*(l+a*T(kk)));
      	SUM_K=SUM_K+xxx'*xxx;
   	end
    delta_2=1/randraw('gamma', [0,1/(HP.chi_tilde_delta_2+1/2*SUM_K),...
     	(HP.xi_tilde_delta_2+1/2*sum(M_k))], [1,1]);    
    clear SUM_K
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(r|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 69-72 in Piecuch_model_description.pdf
    V_R=0; PSI_R=0;
    for kk=1:K
     	if kk==1
         	V_R=V_R+((y_0-b*T0)')*invSigmaMat*(y(:,1)-b*T(1));
         	PSI_R=PSI_R+((y_0-b*T0)')*invSigmaMat*(y_0-b*T0);
        else
         	V_R=V_R+((y(:,kk-1)-b*T(kk-1))')*invSigmaMat*(y(:,kk)-b*T(kk));
          	PSI_R=PSI_R+((y(:,kk-1)-b*T(kk-1))')*invSigmaMat*(y(:,kk-1)-b*T(kk-1));
        end        
   	end
    PSI_R=inv(PSI_R);

    sample=normrnd(PSI_R*V_R,sqrt(PSI_R));
    if sample>HP.v_tilde_r;
        sample=HP.v_tilde_r;
    end
    if sample<HP.u_tilde_r
        sample=HP.u_tilde_r;
    end   
    r=sample;
    clear V_R PSI_R sample

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(sigma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 86-89 in Piecuch_model_description.pdf
    SUM_K=0;
    for kk=1:K
     	if kk==1
         	DYKK=[];
          	DYKK=y(:,1)-r*y_0-(T(1)-r*T0)*b;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        else
         	DYKK=[];
           	DYKK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        end
    end
   	sigma_2=1/randraw('gamma', [0,1/(HP.chi_tilde_sigma_2+1/2*SUM_K),...
     	(HP.xi_tilde_sigma_2+N*K/2)], [1,1]);
   	clear SUM_K DYKK
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(phi|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 73-76 in Piecuch_model_description.pdf
    Phi_now=log(phi);
    Phi_std=0.05;
    Phi_prp=normrnd(Phi_now,Phi_std);
    R_now=C.*exp(-exp(Phi_now)*D);
    R_prp=C.*exp(-exp(Phi_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
        
   	for kk=1:K
      	if kk==1
         	DYYK=y(:,1)-r*y_0-(T(1)-r*T0)*b;
          	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
          	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        else
         	DYYK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b;
         	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
         	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        end
    end
        
 	ins_now=-1/(2*HP.zeta_tilde_phi_2)*(Phi_now-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_phi_2)*(Phi_prp-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-K/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Phi_now=Phi_prp; 
    end
  	phi=exp(Phi_now);
  	clear Phi_now Phi_std Phi_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*

    % redefine matrices since you just updated sigma_2 and phi
    SigmaMat=sigma_2*(C.*exp(-phi*D)); invSigmaMat=inv(SigmaMat);
    SMat=C.*exp(-phi*D); invSMat=inv(SMat);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(l|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 37-40 in Piecuch_model_description.pdf
    V_L=[]; PSI_L=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M,M);
    for kk=1:K
     	SUM_K1=SUM_K1+(selMat(kk).F')*(Z(kk).z-selMat(kk).H*y(:,kk)-...
            selMat(kk).F*a*T(kk));
        SUM_K2=SUM_K2+(selMat(kk).F')*selMat(kk).F;
    end
    V_L=nu/tau_2*ONE_M+1/delta_2*SUM_K1;
    PSI_L=inv(1/tau_2*I_M+1/delta_2*SUM_K2);
    l=mvnrnd(PSI_L*V_L,PSI_L)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(nu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 65-68 in Piecuch_model_description.pdf
    V_NU=[]; PSI_NU=[];
   	V_NU=HP.eta_tilde_nu/HP.zeta_tilde_nu_2+1/tau_2*(ONE_M'*l);
    PSI_NU=inv(1/HP.zeta_tilde_nu_2+M/tau_2);
    nu=normrnd(PSI_NU*V_NU,sqrt(PSI_NU));
    clear V_NU PSI_NU
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tau_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 97-98 in Piecuch_model_description.pdf
    tau_2=1/randraw('gamma', [0,1/(HP.chi_tilde_tau_2+...
     	1/2*(((l-nu*ONE_M)')*(l-nu*ONE_M))),(HP.xi_tilde_tau_2+M/2)], [1,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(gamma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 104-105 in Piecuch_model_description.pdf
    gamma_2=1/randraw('gamma', [0,1/(HP.chi_tilde_gamma_2+...
     	1/2*a'*a),(HP.xi_tilde_gamma_2+M/2)], [1,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(mu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 61-64 in Piecuch_model_description.pdf
    V_MU=[]; PSI_MU=[];
    V_MU=HP.eta_tilde_mu/HP.zeta_tilde_mu_2+ONE_N'*invPiMat*(b+calg*u-calg*wg);
   	PSI_MU=inv(1/HP.zeta_tilde_mu_2+ONE_N'*invPiMat*ONE_N);
    mu=normrnd(PSI_MU*V_MU,sqrt(PSI_MU));
    clear V_MU PSI_MU  
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(pi_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 83-85 in Piecuch_model_description.pdf
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*((b+calg*u-mu*ONE_N-calg*wg)'*invLMat)*(b+calg*u-mu*ONE_N-calg*wg);
    pi_2=1/randraw('gamma', [0,1/(HP.chi_tilde_pi_2+inside2),...
      	(HP.xi_tilde_pi_2+inside1)], [1,1]);
   	clear inside*
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(lambda|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 77-79 in Piecuch_model_description.pdf
    Lambda_now=log(lambda);
    Lambda_std=0.05;
    Lambda_prp=normrnd(Lambda_now,Lambda_std);
    R_now=exp(-exp(Lambda_now)*D);
    R_prp=exp(-exp(Lambda_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
        
    sumk_now=(b+calg*u-mu*ONE_N-calg*wg)'*invR_now*(b+calg*u-mu*ONE_N-calg*wg);
    sumk_prp=(b+calg*u-mu*ONE_N-calg*wg)'*invR_prp*(b+calg*u-mu*ONE_N-calg*wg);
        
 	ins_now=-1/(2*HP.zeta_tilde_lambda_2)*(Lambda_now-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_lambda_2)*(Lambda_prp-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Lambda_now=Lambda_prp; 
    end
  	lambda=exp(Lambda_now);
  	clear Lambda_now Lambda_std Lambda_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*

    % update matrices since you just updated pi_2 and lambda
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat); 
    LMat=exp(-lambda*D); invLMat=inv(LMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(b|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 21-24 in Piecuch_model_description.pdf
   	V_B=[]; PSI_B=[]; SUM_K=ZERO_N;
    for kk=1:K
     	if kk==1
         	SUM_K=SUM_K+(T(1)-r*T0)*(y(:,1)-r*y_0);
        else
          	SUM_K=SUM_K+(T(kk)-r*T(kk-1))*(y(:,kk)-r*y(:,kk-1));
        end
   	end
    V_B=invPiMat*(mu*ONE_N+calg*wg-calg*u)+invSigmaMat*SUM_K;
    PSI_B=inv(invPiMat+invSigmaMat*sum((T-r*[T0 T(1:K-1)]).^2));
    b=mvnrnd(PSI_B*V_B,PSI_B)';
    clear V_B PSI_B SUM_K
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(epsilon_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 99-100 in Piecuch_model_description.pdf
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*(v-u)'*(v-u);
    epsilon_2=1/randraw('gamma', [0,1/(HP.chi_tilde_epsilon_2+inside2),...
      	(HP.xi_tilde_epsilon_2+inside1)], [1,1]);
   	clear inside*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(omega_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 101-103 in Piecuch_model_description.pdf
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*(u-calg*ug-alpha*ONE_N)'*invTMat*(u-calg*ug-alpha*ONE_N);
    omega_2=1/randraw('gamma', [0,1/(HP.chi_tilde_omega_2+inside2),...
      	(HP.xi_tilde_omega_2+inside1)], [1,1]);
   	clear inside*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 80-82 in Piecuch_model_description.pdf
    Rho_now=log(rho);
    Rho_std=0.05;
    Rho_prp=normrnd(Rho_now,Rho_std);
    L_now=exp(-exp(Rho_now)*D);
    L_prp=exp(-exp(Rho_prp)*D);
    invL_now=inv(L_now);
    invL_prp=inv(L_prp);
    sumk_now=0;
    sumk_prp=0;
        
    sumk_now=(u-alpha*ONE_N-calg*ug)'*invL_now*(u-alpha*ONE_N-calg*ug);
    sumk_prp=(u-alpha*ONE_N-calg*ug)'*invL_prp*(u-alpha*ONE_N-calg*ug);
        
 	ins_now=-1/(2*HP.zeta_tilde_rho_2)*(Rho_now-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_rho_2)*(Rho_prp-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_prp;
  	MetFrac=det(L_prp*invL_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Rho_now=Rho_prp; 
    end
  	rho=exp(Rho_now);
  	clear Rho_now Rho_std Rho_prp mat_now mat_prp ins_* sumk MetFrac success_rate L_*

    % redefine matrices because you just updated omega_2 and rho
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat);
    TMat=exp(-rho*D); invTMat=inv(TMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(alpha|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 53-56 in Piecuch_model_description.pdf
    V_ALPHA=[]; PSI_ALPHA=[];
    V_ALPHA=HP.eta_tilde_alpha/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*(u-calg*ug);
   	PSI_ALPHA=inv(1/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*ONE_N);
    alpha=normrnd(PSI_ALPHA*V_ALPHA,sqrt(PSI_ALPHA));
    clear V_ALPHA PSI_ALPHA  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(v|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 33-36 in Piecuch_model_description.pdf
   	V_V=[]; PSI_V=[]; 
    V_V=u/epsilon_2+E'*invDeltaMat*x;
    PSI_V=inv(I_N/epsilon_2+E'*invDeltaMat*E);
    v=mvnrnd(PSI_V*V_V,PSI_V)';
    clear V_V PSI_V
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(u|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 29-32 in Piecuch_model_description.pdf
   	V_U=[]; PSI_U=[]; 
    V_U=calg*invPiMat*(mu*ONE_N+wg-b)+invOmegaMat*(calg*ug+alpha*ONE_N)+v/epsilon_2;
    PSI_U=inv(calg*invPiMat+invOmegaMat+I_N/epsilon_2);
    u=mvnrnd(PSI_U*V_U,PSI_U)';
    clear V_U PSI_U
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(ug|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 45-48 in Piecuch_model_description.pdf
   	V_UG=[]; PSI_UG=[]; 
    V_UG=HP.inv_delta_tilde_ug_2*HP.eta_tilde_ug+calg*invOmegaMat*(u-alpha*ONE_N)+...
        (1/eps_2)*selG'*diagtee'*(diagtee*selG*wg+selD*iota-wye);
    PSI_UG=inv(HP.inv_delta_tilde_ug_2+calg*invOmegaMat+(1/eps_2)*selG'*diagtee'*diagtee*selG);
    ug=mvnrnd(PSI_UG*V_UG,PSI_UG)';
    clear V_UG PSI_UG
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(wg|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 41-44 in Piecuch_model_description.pdf
   	V_WG=[]; PSI_WG=[]; 
    V_WG=HP.inv_delta_tilde_wg_2*HP.eta_tilde_wg+calg*invPiMat*(b+u-mu*ONE_N)+...
        (1/eps_2)*selG'*diagtee'*(wye+diagtee*selG*ug-selD*iota);
    PSI_WG=inv(HP.inv_delta_tilde_wg_2+calg*invPiMat+(1/eps_2)*selG'*diagtee'*diagtee*selG);
    wg=mvnrnd(PSI_WG*V_WG,PSI_WG)';
    clear V_WG PSI_WG

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tee|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 49-52 in Piecuch_model_description.pdf
   	V_TEE=[]; PSI_TEE=[]; 
    V_TEE=(1/eps_2)*diag(selG*(wg-ug))'*(wye-selD*iota)+invXiMat*S;
    PSI_TEE=inv((1/eps_2)*diag(selG*(wg-ug))'*diag(selG*(wg-ug))+invXiMat);
    tee=mvnrnd(PSI_TEE*V_TEE,PSI_TEE)';
    diagtee=diag(tee);
    clear V_TEE PSI_TEE

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(wye|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 13-16 in Piecuch_model_description.pdf
   	V_WYE=[]; PSI_WYE=[]; 
    V_WYE=invGammaMat*H+(1/eps_2)*(diagtee*selG*(wg-ug)+selD*iota);
    PSI_WYE=inv(invGammaMat+(1/eps_2)*I_N_D);
    wye=mvnrnd(PSI_WYE*V_WYE,PSI_WYE)';
    clear V_WYE PSI_WYE

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(beta|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 57-60 in Piecuch_model_description.pdf
    V_BETA=[]; PSI_BETA=[];
    V_BETA=HP.eta_tilde_beta/HP.zeta_tilde_beta_2+(1/kappa_2)*ONE_N_S'*iota;
   	PSI_BETA=inv(1/HP.zeta_tilde_beta_2+N_S/kappa_2);
    beta=normrnd(PSI_BETA*V_BETA,sqrt(PSI_BETA));
    clear V_BETA PSI_BETA  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(iota|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 25-28 in Piecuch_model_description.pdf
    V_IOTA=[]; PSI_IOTA=[];
    V_IOTA=beta/kappa_2*ONE_N_S+(1/eps_2)*selD'*(wye-diagtee*selG*(wg-ug));%nu/tau_2*ONE_M+1/delta_2*SUM_K1;
    PSI_IOTA=inv((1/kappa_2)*I_N_S+(1/eps_2)*selD'*selD);%inv(1/tau_2*I_M+1/delta_2*SUM_K2);
    iota=mvnrnd(PSI_IOTA*V_IOTA,PSI_IOTA)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(kappa_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 93-94 in Piecuch_model_description.pdf
    inside1=[]; inside2=[];
    inside1=N_S/2;
    inside2=1/2*(iota-beta*ONE_N_S)'*(iota-beta*ONE_N_S);
    kappa_2=1/randraw('gamma', [0,1/(HP.chi_tilde_kappa_2+inside2),...
      	(HP.xi_tilde_kappa_2+inside1)], [1,1]);
   	clear inside*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(eps_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see equations 95-96 in Piecuch_model_description.pdf
    inside1=[]; inside2=[];
    inside1=N_D/2;
    inside2=1/2*(wye-diagtee*selG*(wg-ug)-selD*iota)'*(wye-diagtee*selG*(wg-ug)-selD*iota);
    eps_2=1/randraw('gamma', [0,1/(HP.chi_tilde_eps_2+inside2),...
      	(HP.xi_tilde_eps_2+inside1)], [1,1]);
   	clear inside*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now update arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_all_arrays

end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete the burn-in period values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete_burn_in

%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%
save(['bayes_model_solutions/experiment_',save_name,'.mat'],...
    'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
    'PHI','LAMBDA','A','B','ELL','R','Y_0','Y', 'U','UG','WG','V','ALPHA',...
    'RHO','OMEGA_2','EPSILON_2','HP','*DATA','*LON','*LAT',...
    '*NAME','*ERROR','N','K','M','L','D','nn','HOL_*',...
    'KAPPA_2','EPS_2','BETA','IOTA','WYE','TEE','inp')