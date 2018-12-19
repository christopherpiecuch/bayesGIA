%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [V_GPS,S_GPS,X_GPS,Y_GPS,N_GPS,T_GPS] = prepare_gps_data(la1,la2,lo1,lo2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads GPS vertical velocities and standard errors from the ULR6a dataset
% from stations found within the geographic study domain
% INPUT: 
%   la1         Southern latitudinal bounds of study region
%   la2         Northern latitudinal bounds "
%   lo1         Western longitudinal bounds "
%   lo2         Eastern longitudinal bounds "
% OUTPUT:
%   V_GPS       Vertical velocities
%   S_GPS       Standard errors
%   X_GPS       Longitudes
%   Y_GPS       Latitudes
%   N_GPS       Names of GPS sites
%   T_GPS       Length of time series from which V and S were computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V_GPS,S_GPS,X_GPS,Y_GPS,N_GPS,T_GPS] = prepare_gps_data(la1,la2,lo1,lo2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bring in ULR6a GPS VLM trends (and errors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('GPS_ULR6a/ulr6a_vertical_velocities-20171024.mat','Site','Lon','Lat','V_GPS','S_GPS','T_GPS')

% Delete locations not within the box specified above
i = find(~(Lon>lo1&Lon<lo2&Lat>la1&Lat<la2));
Lat(i)=[];
Lon(i)=[];
S_GPS(i)=[];
Site(i,:)=[];
V_GPS(i)=[];
T_GPS(i)=[];

% in case you're doing the east coast, delete pesky open-ocean sites
for n=numel(Lat):-1:1
   if strcmp(Site(n,:),'EXU0')
       Lat(n)=[];
       Lon(n)=[];
       Site(n,:)=[];
       V_GPS(n)=[];
       S_GPS(n)=[];
       T_GPS(n)=[];
   end
   if strcmp(Site(n,:),'SCUB')
       Lat(n)=[];
       Lon(n)=[];
       Site(n,:)=[];
       V_GPS(n)=[];
       S_GPS(n)=[];
       T_GPS(n)=[];
   end
   if strcmp(Site(n,:),'NAS0')
       Lat(n)=[];
       Lon(n)=[];
       Site(n,:)=[];
       V_GPS(n)=[];
       S_GPS(n)=[];
       T_GPS(n)=[];
   end
   if strcmp(Site(n,:),'BRMU')
       Lat(n)=[];
       Lon(n)=[];
       Site(n,:)=[];
       V_GPS(n)=[];
       S_GPS(n)=[];
       T_GPS(n)=[];
   end
end

% One issue with the GPS data is that, in addition to their sparseness,
% etc., the sites can be very close to one another .. or literally right on
% top of one another. This doesn't bode well for all the matrix inverses we
% have to compute in the main code (i.e., they don't exist). To address
% this, I employ a two step algorithm. First, I concatonate stations that
% are literally right on top of one another. Second, I do a subsequent
% search for, given the remainder, are there any pairs that are "really
% close" (defined by some critical distance). If there are, I only keep the
% longer of the records.

% see if there are any with the *EXACT* same location
% if so, combine assuming independent errors
for n=numel(Lat):-1:1
   i=[]; i=find(Lat==Lat(n)&Lon==Lon(n));
   i(i==n)=[];
   if ~isempty(i)
       Lat(n)=[];
       Lon(n)=[];
       Site(n,:)=[];
       V_GPS(i)=0.5*V_GPS(i)+0.5*V_GPS(n);
       S_GPS(i)=sqrt(1/4*S_GPS(i)^2+1/4*S_GPS(n)^2);
       T_GPS(i)=T_GPS(i)+T_GPS(n);
       V_GPS(n)=[];
       S_GPS(n)=[];
       T_GPS(n)=[];
   end
end

X_GPS = Lon;
Y_GPS = Lat;
N_GPS = Site; % keep the names
V_GPS = 1e-3*V_GPS; % convert to m/yr
S_GPS = 1e-3*S_GPS; % convert to m/yr

return
