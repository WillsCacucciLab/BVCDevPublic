function [M] = makeCanonicalModelMaps_v2(varargin)
% Make the complete set of BVC model maps (for given phi and d ranges) for
% 'canonical' environment shapes, as follows:
%
% 'hp' - Square environment 25x25 pixels
% 'barrier N-S'
% 'barrier W-E'
%
% etc .. FINISH ME!
%
% Note that all the env shape dimensions are hard-coded in this function.
% 06/20 LM: changed the barrStr for all barr trials as format wasn't consistent with that of baseline 
% baseline: xl,yl, xl,yu, 0; xl,yu, xu,yu, 0; xu,yu, xu,yl, 0; xu,yl, xl,yl, 0;
% barr before change: xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1; xl,yl, xl,yu, 1;


prms.d_list                   = 1:16;                % Distance tunings. Defined in map bins, not cm!
                                                     % .. note that 1 is minimum meaningful value, as edge of env defined 1 bin from wall.
                                                     % .. Max of d_list is max possible over all envs, in smaller envs where d>(mapSize/2),
                                                     % .. model will be all NaN.
prms.phi_list                 = deg2rad( 0:6:354 );  % Angular tunings. In radians.
prms.binSize                  = 2.5;                 % Bin size in cm. Needed only to translate the radial tuning parameteters from Hartley et al 2000 (defined in cm) into map bin units.
prms.kappa                    = 25;      % \
prms.beta                     = 183;     %  - This is the tuning parameter set. 
prms.sigZeros                 = [6.2   12.2   20.2   30.2];    % (2.4928 : 1 : 6.4928) .^ 2; % sqrt(sig0) varies linearly and sqrt(12.2)=3.4928

prms.save                     = 1; % y/n
% Allow caller input for params %
if ~isempty(varargin)
    for ii=1:2:(length(varargin)-1)
        prms.(varargin{ii}) = varargin{ii+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1). Make the coordinate sets that define the shapes of the canonical envs.
%      The important output of this this section is the cell array 'canonBarrMaps', the logical arrays
%      that define the shape of each environment. (Along with text key, 'canonicalEnvs').
n = 0;
% 1a. Baseline map (open field square) %
n = n + 1;
mapDefs(n).env      = 'hp';
mapDefs(n).mapSize  = [27 27];
mapDefs(n).barrStr  = [1,1, 1,27, 0; 1,27, 27,27, 0; 27,27, 27,1, 0; 27,1, 1,1, 0];
mapDefs(n).barrCirc = [];

% % 1b. Barrier N=S map %
n = n + 1;
mapDefs(n).env      = 'barrier N-S';
mapDefs(n).mapSize  = [27 27];
xl=13.5;  xu=14.49;  yl=4;  yu=24;
mapDefs(n).barrStr  = [mapDefs(1).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1];
mapDefs(n).barrCirc = [];
% barrMap(4:24, 14) = 1;

% 1c. Barrier W-E map %
n = n + 1;
mapDefs(n).env      = 'barrier W-E';
mapDefs(n).mapSize  = [27 27];
xl=4;  xu=24;  yl=13.5;  yu=14.49;
mapDefs(n).barrStr  = [mapDefs(1).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1];
mapDefs(n).barrCirc = [];
% barrMap(14, 4:24) = 1;

% 1d. Barrier N=S map (short fat version) %
n = n + 1;
mapDefs(n).env      = 'barrier N-S short';
mapDefs(n).mapSize  = [27 27];
xl=13;  xu=15;  yl=6;  yu=22;
mapDefs(n).barrStr  = [mapDefs(1).barrStr;xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1 ];
mapDefs(n).barrCirc = [];


% 1e. Barrier W-E map (short fat version) %
n = n + 1;
mapDefs(n).env      = 'barrier W-E short';
mapDefs(n).mapSize  = [27 27];
xl=6;  xu=22;  yl=13;  yu=15;
mapDefs(n).barrStr  = [mapDefs(1).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1 ];
mapDefs(n).barrCirc = [];

% 1f. Circle (aka CCE). This has been scaled to be exactly 32 pixels diameter, so sits in the centre of a 34x34 array. %
n = n + 1;
mapDefs(n).env      = 'CCE';
mapDefs(n).mapSize  = [34 34];
mapDefs(n).barrStr  = [];
mapDefs(n).barrCirc = [17.5 17.5 16.5 0];

% 1g. Baseline map (open field square) for Moser Data, so 70cm x 70cm %
n = n + 1;
mapDefs(n).env      = 'hp70';
mapDefs(n).mapSize  = [30 30];  % 70cm square => 28 bins @2.5cm/bin
mapDefs(n).barrStr  = [1,1, 1,30, 0; 1,30, 30,30, 0; 30,30, 30,1, 0; 30,1, 1,1, 0];
mapDefs(n).barrCirc = [];
hp70Ind             = n;  % Keep a record of this, for making Moser barrier maps.

% the Moser barrier measures 35x1cm according to Bjerknes et al.; so in
% bins: 0.4x14 - however we need to make the barrier twice the size as the
% environment has even number of bins so the barrier could be on either
% side of the centre
% 1g. Barrier VT map (Moser Data, vertical barrier abutting top wall) %
n = n + 1;
mapDefs(n).env      = 'barrier VB'; % LM bugfix - wrong definition here: mapDefs(n).env = 'barrier VT';
mapDefs(n).mapSize  = [30 30];  
xl=14.8;  xu=16.2;  yl=16;  yu=30; % xl=14;  xu=17;  yl=15;  yu=30;
mapDefs(n).barrStr  = [mapDefs(hp70Ind).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1 ];
mapDefs(n).barrCirc = [];

% 1h. Barrier VB map (Moser Data, vertical barrier abutting bottom wall) %
n = n + 1;
mapDefs(n).env      = 'barrier VT'; % LM bugfix - wrong definition here: mapDefs(n).env = 'barrier VB';
mapDefs(n).mapSize  = [30 30];  
xl=14.8;  xu=16.2;  yl=1;  yu=15; % xl=14;  xu=17;  yl=1;  yu=15;
mapDefs(n).barrStr  = [mapDefs(hp70Ind).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1];
mapDefs(n).barrCirc = [];

% 1i. Barrier HR map (Moser Data, horizontal barrier abutting right wall) %
n = n + 1;
mapDefs(n).env      = 'barrier HR';
mapDefs(n).mapSize  = [30 30];  
xl=16;  xu=30;  yl=14.8;  yu=16.2; % xl=15;  xu=30;  yl=14;  yu=17;
mapDefs(n).barrStr  = [mapDefs(hp70Ind).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1];
mapDefs(n).barrCirc = [];

% 1j. Barrier HL map (Moser Data, horizontal barrier abutting left wall) %
n = n + 1;
mapDefs(n).env      = 'barrier HL';
mapDefs(n).mapSize  = [30 30];  
xl=1;  xu=15;  yl=14.8;  yu=16.2; % xl=1;  xu=15;  yl=14;  yu=17;
mapDefs(n).barrStr  = [mapDefs(hp70Ind).barrStr; xl,yl, xl,yu, 1; xl,yu, xu,yu, 1; xu,yu, xu,yl, 1; xu,yl, xl,yl, 1 ];
mapDefs(n).barrCirc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2). Make the complete sets of BVC models. 
[d_grid, phi_grid, sz_grid]  = meshgrid(prms.d_list, prms.phi_list, prms.sigZeros);    
modelBVCs                    = cell(1, length(mapDefs));  % Pre-allocate storage for these BVCs maps    
parfor ii=1:length(mapDefs) % Loop is PARFOR compatible (may have been switched to FOR for debugging). 
    
    % Generate the BVC maps for this trial: run through all d and phi combinations.
    modelBVCsTemp = nan( [mapDefs(ii).mapSize, numel(d_grid)] );
    for jj = 1:numel(d_grid)  % jj=iterator over all possible models (d, phi and sigma zero)
        
        % First check if proposed d is greater than env side length/2, for this env.
        % If so, don't make model map, leave an entry but set as NaN.
        if d_grid(jj)>ceil((mapDefs(ii).mapSize-2)/2)
            modelBVCsTemp(:,:,jj) = nan( mapDefs(ii).mapSize );
            continue
        end
        
        % If above not a problem, make model map.
        modelBVCsTemp(:,:,jj) = makeBVCMap_v2( d_grid(jj), phi_grid(jj), mapDefs(ii).mapSize, mapDefs(ii).barrStr, mapDefs(ii).barrCirc, prms.binSize, [prms.beta, sz_grid(jj), prms.kappa] );
    
    end
    modelBVCs{ii} = modelBVCsTemp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign output %
for jj = 1:numel(d_grid)
    [M.phiKey(jj), M.dKey(jj), M.szKey(jj)] = ind2sub( size(d_grid), jj );   % 'Key' for each map will be an index, refernce to 'd_list' and 'phi_list' to get actual value. 
end
M.modelMaps = modelBVCs;
M.mapDefs   = mapDefs;
M.d_list    = prms.d_list;
M.phi_list  = prms.phi_list;
M.sz_list   = prms.sigZeros;
M.tuning    = [prms.beta, 1, prms.kappa];  % '1' here is a dummy to keep format, used to be sigZero

if prms.save
    save( ['CanBVCmaps_' date '.mat'], 'M' );
end
end






