%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ClusMult = determine_clusters(Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch et al., 2018, Origin of spatial variation in United
% States East Coast sea level trends during 1900-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is specific to the spatial covariance structure of sea level
% variation along the U.S. East Coast.  Given a set of N locations, this
% function creates an [NxN] matrix ("C") filled with zeros and ones.  The
% values C_ij equal 1 iff locations i and j are both *either* north or
% south of Cape Hatteras, but 0 if locations i and j are on opposite sites
% of hatteras.
% INPUT: 
%   Y           Latitudes of locations
% OUTPUT:
%   ClusMult    The "C" matrix defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ClusMult = determine_clusters(Y)

% This is a script where you can, in addition to the actual tide gauge
% locations, add in any other geographic sites (target locations) where you
% want the model to impute sea level values. Right now I've written this
% code more or less as a place holder, inserting just one extra point, so
% to be sure the code works and estimates the process at sites with no data
% at all. However, for more realistic applications, you'd enter many more
% points, e.g., a regular grid along the coast, e.g., for comparison with
% climate models or computing "true" spatial averages along the coast.

% This is a script where you will flag which cluster a given site falls
% within. Right now this is more or less a place-holder script, designed
% for application to the North American East Coast and the "divide" at Cape
% Hatteras. Thus, sites are in one cluster or the other based on their
% latitude relative to Hatteras. For more general application, you will
% need to expand and generalize this code. But for now, it's just here so
% that I can make sure the larger body of the code works with clustered
% regions.

Y_HATTERAS = 35.25;
BLOCK(find(Y<Y_HATTERAS)) = 1;
BLOCK(find(Y>Y_HATTERAS)) = 2;

ClusMult = BLOCK'*BLOCK;
ClusMult(ClusMult==2)=0;
ClusMult(ClusMult~=0)=1;

return