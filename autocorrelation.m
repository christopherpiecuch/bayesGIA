%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rk=autocorrelation(Y,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample autocorrelation of a time series Y(t) at lag k
% INPUT: 
%   Y           Time series
%   k           lag
% OUTPUT:
%   rk          k-lag sample autocorrelation of Y(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rk=autocorrelation(Y,k)

n=numel(Y);
kp1=k+1;
Ybar=mean(Y);

top=0; bot=0;

for t=kp1:n
    top=top+(Y(t)-Ybar)*(Y(t-k)-Ybar);
end
for t=1:n
    bot=bot+(Y(t)-Ybar)^2;
end

rk=top/bot;

return