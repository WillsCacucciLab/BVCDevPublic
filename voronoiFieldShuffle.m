function [mapsShuf] = voronoiFieldShuffle(map_unsm,map_pos,nRnd,cellID,trID)
% Randomise data by voronoi-segmenting rate map and randomly rearranging segments.
%
%       [mapsShuf] = voronoiFieldShuffle(map_unsm, map_pos, n_shuffle)
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Map/field segmentation - this only needs to happend once per cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map_spk  = round( map_unsm .* map_pos );

% Smooth before peak finding:
map      = rates_adaptivesmooth( map_pos, map_spk, 200 );

% Remove NaNs from map.
% map = inpaintn(map); TW removed this as part of adapting for circles

% calculate a sort of 'prominence' measure: for every bin, the rate ratio between
% bin plus all 8-connected bins, and the surrounding 18 bins that would make up a 5x5 square.
[cen, surr]     = deal(ones(5,5));
cen([1 5],:)    = 0;   cen(:, [1 5]) = 0;
surr(2:4,2:4)   = 0;
visMask         = double( ~isnan( map_pos ) );
map(isnan(map)) = 0;   % Added this for circles - can't have nans in map array - previously interpolated them with rate values in square
cenConv         = conv2( map, cen ) ./ conv2( visMask, cen );
surrConv        = conv2( map, surr ) ./ conv2( visMask, surr );
prom            = cenConv ./ surrConv;
prom            = prom( (1:size(map,1))+2, (1:size(map,1))+2 );

% Find local maximae %
lMaxMask            = imregionalmax( map );
lMaxInd             = find( lMaxMask );
[lMaxY,lMaxX]       = ind2sub( size(map), lMaxInd );


% Sort by prominence measure, and take best N (8?)
promPks       = prom(lMaxMask);
[~,promOrder] = sort(promPks,'descend');
nPks          = min( [length(lMaxX), 8] );
lMaxY         = lMaxY( promOrder(1:nPks) );
lMaxX         = lMaxX( promOrder(1:nPks) );


% Pad the peak list with a bounding box, to avoid having voronoi regions
% bounded at infinity. Any regions which *are* bounded at infinity will be
% part of the bounding box, and can be ignored. Need to place bounding box
% far from map, otherwise get 'edge effects' - voronoi regions change.
lb    = -100;
ub    = 125;
lMaxY = [lMaxY; lb; lb; ub; ub];
lMaxX = [lMaxX; lb; ub; lb; ub];


% Voronoi via DT function gets vertices and cells:
DT    = delaunayTriangulation(lMaxX,lMaxY);
[V,R] = voronoiDiagram(DT);

% Get zones from map:
nFdZn = 0;
for itZn=1:length(R)
    if any(R{itZn}==1);  continue;  end
    nFdZn         = nFdZn+1;
    znMask{nFdZn} = poly2mask( V( R{itZn}, 1 ), V( R{itZn}, 2 ), size(map,1), size(map,2) ); %%%% TODO - hard coded map size here - how to fix for circles?
end


% Prep the zones for random shifts: get the actual rate values for the zone, 
% crop tight, sort by mean firing rate:
znMR       = nan(nFdZn,1);
znRates    = repmat( {nan(size(map))}, nFdZn, 1);
for itZn = 1:nFdZn
    znRates{itZn}( znMask{itZn} ) = map_unsm( znMask{itZn} );
    znRates{itZn}( isnan(map_pos) )   = nan;  % Inserted for circles - explicitly set non visited rates to nan
    znMR(itZn)                    = mean( znRates{itZn}, 'all', 'omitnan');
    [znR,znC]                     = ind2sub( size(map), find(znMask{itZn}));
    znRates{itZn}                 = znRates{itZn}( min(znR):max(znR), min(znC):max(znC) );
end
[~,rateOrder] = sort(znMR,'descend');
znRates       = znRates( rateOrder );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the randomised bit, to be repeated many times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prep for randomising - pre-rotate all zones to all possible rotations.
% (saves time if N repeats > N poss rotations).
rotList = 0:6:354;
rotZns  = cell(length(rotList),length(znMask));
for itRt=1:length(rotList)
    for itZn=1:length(znMask)
        rotZns{itRt,itZn}             = imrotate( znRates{itZn}, rotList(itRt) );  
        rotMask                       = imrotate( ones(size(znRates{itZn})), rotList(itRt) );  % Need to set padding created by
        rotZns{itRt,itZn}(rotMask==0) = nan;                                                   % rotation to Nan (func makes 0).
    end
end
[nMapR,nMapC] = size(map);
nPossRot      = length(rotList);
isVisited     = map_pos > 0;
mapsShuf      = nan(nMapR,nMapC,nRnd);  % Pre-assign output.

% Reset the random number generator for each map, to create repeatable output
rng('default')

% Random shift and rotate
for itRd = 1:nRnd
    
    try
        
        rndMap       = nan(nMapR,nMapC);
        overflowBins = [];
        for itZn = 1:length(znMask)
            
            % Get rotated zone.
            znRot                           = rotZns{  randi(nPossRot) , itZn };
            
            % Place zone within random map: if the rotated zone is larger than the map in either dimension
            % (an unusual case, for long zones which get rotated), then crop to map size in that dimension
            % before proceeding (crop at a random point). Any cropped bins containing actual rate data are
            % put into the 'overflow' list already (to be assigned bin-by-bin to spare slots, after zone
            % placement.
            if size(znRot,1) > nMapR
                logCrInd                                               = false( size(znRot,1), 1 );
                logCrInd(  (1:nMapR)+randi( size(znRot,1)-nMapR )-1  ) = true;
                overflowBins                                           = [overflowBins; reshape( znRot(  ~logCrInd,  :   ), [], 1 )];
                znRot                                                  = znRot( logCrInd,  :   );
            end
            if size(znRot,2) > nMapC
                logCrInd                                               = false( 1, size(znRot,2) );
                logCrInd(  (1:nMapC)+randi( size(znRot,2)-nMapC )-1  ) = true;
                overflowBins                                           = [overflowBins; reshape( znRot(  :,   ~logCrInd  ), [], 1)];
                znRot                                                  = znRot(   :,   logCrInd  );
            end
            % After this, place the rotated zone in a random place, but where it will fit in the whole map.
            rndTL_R  = randi( nMapR - size( znRot, 1) + 1 ) - 1;   % +1 then -1 because RANDI input must be >=1
            rndTL_C  = randi( nMapC - size( znRot, 2) + 1 ) - 1;   %    "
            tmpZnMap = nan(nMapR,nMapC);
            tmpZnMap( (1:size(znRot,1))+rndTL_R, (1:size(znRot,2))+rndTL_C ) = znRot;
            
            % Put together this rotated, shifted zone with the other zones, checking for overlaps.
            % Overlapping pixels get put on the overflow list, to be assigned (bin-by-bin) later.
            % Note that zones earlier in the list get placed first and have priority to remain 'coherent'.
            % As zones are sorted by mean rate, should mean that firing fields have more chance of remaining
            % coherent, background firing more likely to go in overflowBins.
            % In the two lines below, 'isnan(rndMap)' is testing whether rate data has already been assigned 
            % to a position, following placement of other zones. 'isVisited' is a test of whether bins were occupied
            % in the actual data.
            goodToPlaceMsk         = isnan( rndMap ) & isVisited & ~isnan(tmpZnMap);
            overflowMsk            = (~isnan( rndMap ) | ~isVisited) & ~isnan(tmpZnMap);
            rndMap(goodToPlaceMsk) = tmpZnMap(goodToPlaceMsk);
            overflowBins           = [overflowBins; tmpZnMap(overflowMsk)];
            
        end
        
        % Now need to place the overflow bins. There is almost always a small mismatch between N
        % overflow bins and free places in rndMap, as after rotation the zones are not guaranteed
        % to have the same number of non-nan bins. If too many overflow bins, choose a random subset.
        % If not enough, choose a random subset to double up.
        overflowBins = overflowBins( ~isnan(overflowBins) );  % There will be some NaN bins from when too-big rotated zones get cropped.
        freeBinsInd  = isnan(rndMap) & isVisited;
        nPlaces      = sum(  freeBinsInd  ,    'all'    );
        nBins        = length( overflowBins );
        if nPlaces > nBins
            repInd       = randperm(nBins, nPlaces-nBins);
            overflowBins = [overflowBins; overflowBins(repInd)];
        elseif nPlaces < nBins
            overflowBins = overflowBins( randperm(nBins, nPlaces) );
        end
        rndMap(freeBinsInd) = overflowBins;
        
        % Smooth random map
        rndMap(~isVisited) = 0;
        filtMap            = imfilter( rndMap, ones(5) );
        filtVis            = imfilter( visMask, ones(5) );
        rndMap             = filtMap ./ filtVis;
        rndMap(~isVisited) = nan;
        
        % Assign to output
        mapsShuf(:,:,itRd) = rndMap;
    
    catch ME
        
        errorLog.exception = ME;
        errorLog.cellID    = cellID;
        errorLog.trID      = trID;
        errorLog.itRand    = itRd;
        save( ['vs_errorLog_' cellID '_tr' num2str(trID) '.mat'], 'errorLog' );
        
        
    end

    
end

