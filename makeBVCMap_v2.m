function [bvcMap] = makeBVCMap_v2( dist, phi, mapSize, barrStr, barrCirc, binSize, tuning, varargin )
% Make BVC map, v2: barriers specified as vectors, represented at high spatial sampling frequency.
%
% Two ways to use. (1) Model map way (this is the normal way): 
%
%           bvcMap = makeBVCMap_v2( dist, phi, mapSize, barrStr, barrCirc, binSize, tuning )
%
% barrStr format:  [x1 y1, x2 y2, barrID;   etc   ] for each straight line wall segment.
%                     
% barrCirc format: [xCen, yCen, radius, barrID;   etc  ] for each circle.
%
% The units used for: dist, mapSize, barrStr and barrCirc are *map bins*. The units of argument 'binSize' 
% are cm: this is used only to define the tuning parameters in 'tuning' in terms of bins. (Tuning parameters
% are defined following Hartley et al, and units are cm).
%
% If tuning=[], tuning parameters will default to Hartley et al values. 
%
% 'barrID' is 0 for all walls forming the outer wall of the arena, otherwise 1, 2 etc for each 
%  separate physical barrier within the arena (including if it is abutting an outer wall).
%
% Usage (2) - simulate a BVC using real path data and a given approximate peak rate (firing is Poisson).
%
%       bvcMap = makeBVCMap_v2( dist, phi, mapSize, barrStr, barrCirc, binSize, tuning, XY, meanRate )
%
% Generally as before, but now supply X and Y, which are actual path data, and meanRate, which is the 
% desired approximate mean rate of the BVC. By-sample spiking is random poisson. <XY format is (nPosSamp,x:y)>.
%
% You need to make the units of X,Y and dist, mapSize, barrStr and barrCirc match up, still. If you are providing 
% X,Y at camera pixel resolution, for example, then you enviroment shape definitions need to be in the same units.

 
% Written as a replacement for Ben Towse's functions 'apply_BV_model' and 'predictBVFiring'. The problem with those functions were inaccuracies in the
% radial sampling of the BVC function (basically, the problem stemmed from evaluating the function at bin centres, as you came right up against a wall, 
% the angular distance between bin centres was so large that you fail to sample properly the BVC tuning curve. The tuning curve 'fit inbetween' the points
% at which it was being sampled).
% This function gets around that by sampling the position of the wall at a high spatial rate (10x bin spatial frequency), but at the price of a somewhat
% clunky methods of defining wall positions (need to give coords of each straight line) and occlusion detection (need to make sure BVC function is only 
% evaluated at closest wall, not occluded ones).
%
% Tom Wills, Jan 2018.

if isempty(tuning)
    beta    = 183 / binSize;
    sigZero = 12.2 / binSize;
    kappa   = 25;
else
    beta    = tuning(1) / binSize;
    sigZero = tuning(2) / binSize;
    kappa   = tuning(3);
end
sigRad    = sigZero*(  (dist/beta)+1  );

if isempty(varargin)
    simMode = 0;
else
    simMode = 1;
    actX    = varargin{1}(:,1);   actY  = varargin{1}(:,2);     simMeanRate  = varargin{2};
    nanInd  = isnan(actX) | isnan(actY);
    actX    = floor( actX(~nanInd) );   % X and Y should already be min = 1, so round is 
    actY    = floor( actY(~nanInd) );   % correct here.
%     actX    = actX( 1:1500 );
%     actY    = actY( 1:1500 );
end

wallSegOverSamp = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Make the wall and barrier XY coord list. (convert corner coords or centre/radius specs to lines).
%     The main output from this section is 'barrCoords', with dimensions ( nWallPoints, 4 ), where the format of each line is:
%
%           [x, y, barrID, wallID]
%
%     'barrID', is supplied by user, corresponds to individual physical barrier or object. barrID=0 is 
%     reserved for the outer wall.
%     'wallID', is 0 for all outer walls, and then increments 1 per actual wall segment otherwise.
%     Used for removing non-visible points. (Note that this system assumes that outer walls never block view of themselves, though).
barrCoords  = ones( 0, 4 );
wallIDCount = 1;
% 1a. Straight wall segments
for ii=1:size( barrStr,1 )
    barrLength = sqrt( diff( barrStr(ii,[1 3]) ).^2 + diff( barrStr(ii,[2 4]) ).^2  );  % First work out length then assign number of points of line,
    nPoints    = round( barrLength*wallSegOverSamp );                                                % this way, keep even 1/10 grid spacing, even for diagonals.
    barrXTemp  = linspace( barrStr(ii,1), barrStr(ii,3), nPoints )';
    barrYTemp  = linspace( barrStr(ii,2), barrStr(ii,4), nPoints )';
    if barrStr(ii,5)==0
        wallIDTemp = zeros(size(barrXTemp));
    else
        wallIDTemp = ones(size(barrXTemp)) .* wallIDCount;    wallIDCount = wallIDCount + 1;
    end
    barrCoords = cat(1, barrCoords, [barrXTemp, barrYTemp, ones(size(barrXTemp)).*barrStr(ii,5), wallIDTemp]);
end
% 1b. Circle wall segments.
%     TODO - need to control that user input is OK, that there is an outer wall (or it fails gracefully if not), that if there are barrStr and barrCirc
%     inputs, the barrIDs are not conflicting or overlapping.
for ii=1:size( barrCirc, 1 )
    nPoints = round( 2*pi*barrCirc(ii,3)*wallSegOverSamp );  % nPoints is circumference (in bins) x10.
    [barrXTemp, barrYTemp] = pol2cart( linspace( 2*pi/nPoints, 2*pi, nPoints )', ones(nPoints,1).*barrCirc(ii,3) );  % First angular binin circle is not 0, so 0/2pi are not double counted.
    barrXTemp              = barrXTemp + barrCirc( ii, 1 );
    barrYTemp              = barrYTemp + barrCirc( ii, 2 );
    if barrCirc(ii,4)==0
        wallIDTemp = zeros(size(barrXTemp));
    else
        wallIDTemp = ones(size(barrXTemp)) .* wallIDCount;    wallIDCount = wallIDCount + 1;
    end
    barrCoords = cat(1, barrCoords, [barrXTemp, barrYTemp, ones(size(barrXTemp)).*barrCirc(ii,4), wallIDTemp]);
end
nBarrWalls = wallIDCount - 1;
% 1c. If wall segments are continuous, there will be an overlapping point where they meet: remove these by getting uniqueXY rows.
[barrXYUni, ind] = unique( barrCoords(:,1:2), 'rows', 'stable' );  % doinf unique on barrCoords without a row index only gets the corners withon a wall ID, i.e. for the outer wall.
barrCoords       = [barrXYUni, barrCoords(ind,3:4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Some more prep work on the map: (a) work out which map bins the barriers/walls fall in, 
%     so as we can exclude these bins from the analysis, (b) find the convex hull of the outer 
%     wall (that with barrID assigned to 0), and exclude bins outside of this from the analysis.
if ~simMode
    visMap = true( mapSize );   % 'visMap' will be a logical mask of the bins in the map for which we will create a BVC response.
else
    visMap                = false( mapSize );
    [visInd,~,pos2VisInd] = unique(  sub2ind(mapSize,actY+1,actX+1)  );
    visMap(visInd)        = true;
end
% 2a. Find barrier bins in map
barrBinsInd          = sub2ind( mapSize, round(barrCoords(:,2)), round(barrCoords(:,1)) );  % Bin<->coord convention is that integers represent bin centres, so 'round' is the correct operation. 
visMap( barrBinsInd) = false;
% 2b. Exclude bins outside outer wall (a procedure which should be avoidable for square, but is essential for circle), 
%     and also those that are inside barriers of a >1bin thickness.
% To save time, don't do this for squares with no internal barriers. As long as user input is 'sensible'
% (walls run along outer bin of map), then this bit isn't necessary.
barrIDList = unique( barrCoords(:,3) );
if length(barrIDList)>1   ||  ~isempty( barrCirc )  
    [X,Y]      = meshgrid( 1:mapSize(2), 1:mapSize(1) );
    for ii=1:length( barrIDList )
        indBarr    = barrCoords(:,3) == barrIDList(ii);
        coordsTemp = barrCoords( indBarr, : );
        indCH   = convhull( coordsTemp(:,1), coordsTemp(:,2) );
        inMask  = inpolygon( X, Y, coordsTemp(indCH,1), coordsTemp(indCH,2) );
        if barrIDList(ii)==0
            visMap(~inMask) = false;  % If the barrier in question is the outer wall, then set the visMap *outside* to false.
        else
            visMap(inMask) = false;   % If it is an internal barrier, then set the visMap *inside* the barrier to false.
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Now make the BVC map: loop through every visited bin (as defined in 'visMap'), calculate 
%     the distance and direction to all ofthe visible wall coordinates, and evaluate the BVC
%     function at these points.
[visY,  visX]  = ind2sub( mapSize, find(visMap) );
bvcMap         = nan( mapSize );
for ii=1:length(visX)  % ii=iterator over visited bins

    % 3a. Get a circularly sorted list of barrier position coords (polar), *relative to visited bin ii* %
    [Th, R]           = cart2pol( barrCoords(:,1)-visX(ii), barrCoords(:,2)-visY(ii) );
    [ThSort, sortInd] = sort( Th, 'ascend' );
    RSort             = R( sortInd );
    wallIDSort        = barrCoords( sortInd, 4 );
    
    % 3b. If we are looking stright down a wall from the end (i.e. wall is self occluding), only take the closest segment.
    % Code here assumes that there is only one self-occluding wall at any one time, otherwise won't work.
    diffTh                = diff( ThSort );
    if sum( diffTh==0 ) > 10      % There might always be a few pairs of identical thetas (from different walls), but continuous self-occluding wall must have more than 10. 
        constThInd         = diffTh==0 & diff(wallIDSort)==0;  % The second test here is to make sure it is continuous theta from a continuous wall.
        RSortTemp          = RSort;
        RSortTemp( ~constThInd ) = nan;
        [~, closestSegInd] = nanmin( RSortTemp );
        constThInd(  closestSegInd  ) = false;
        ThSort         = ThSort( ~constThInd );
        RSort          = RSort( ~constThInd );
        wallIDSort         = wallIDSort( ~constThInd );
    end
    
    % 3c. When there is a barrier in the middle of the box, need a routine to get the *closest* wall at each angle %
    %     The approch is based on the assumption that walls forwhich the closest point to the rat is the least occulde
    %     others which have the closest point further away. Ths is fine for a simple rectangle barrier, but will not
    %     work for many situations with more complex shapes, or with two barriers.
    % First, get the closest distance to the rat for all of the internal (barrier) walls. %
    closestWallDists = nan(1,nBarrWalls);
    for jj=1:nBarrWalls
        closestWallDists(jj) = min( RSort( wallIDSort==jj ) );
    end
    [closestWallDistsSorted,closestWallOrder]      = sort(closestWallDists,'ascend');
    % Debug stop %
    if 0 && visY(ii)==3 && visX(ii)==14
        disp('dum');
    end
    % Then, starting from the closest wall, remove points from all other walls that lie within the extent of that wall %
    for jj=1:length(closestWallOrder)
        wallNumInd   = find( wallIDSort==closestWallOrder(jj) );
        if length(wallNumInd) <= 1
            continue  % In the case that wall N has already been completely removed, or having just one segment left (which cant occlude anything).
        end        
        % The only tricky bit here comes from deciding if the wall is lying across the angular wrap point (as the wall segments
        % will never be a completely coherent block, they will always be interspersed with segments of other occluding or occluding wall). 
        wallGaps      = diff( wallNumInd );
        wallExtentInd = false( size(wallIDSort) );
        gapTol        = 300;
        if max(wallGaps)<gapTol
            wallExtentInd( wallNumInd(1) : wallNumInd(end) ) = true;    % In the case of a coherent chunk in the middle of the circular FOV ..
        else
            wallExtentInd( 1 : wallNumInd(    find(wallGaps>gapTol,1,'first')    ) )  = true;  % If the wall lies across the wrap point
            wallExtentInd( wallNumInd(   find(wallGaps>gapTol,1,'first')+1   ) : end) = true;
        end
        % Now remove wall segments not from wall N, lying within the extent of wall N %
        segRemoveInd = wallExtentInd & wallIDSort~=closestWallOrder(jj);
        % Check if the closest point is on a corner of two walls: if it is, don't remove points from the other wall that makes up the corner, 
        % so as to preserve the corner (otherwise, if the walls have sorted the 'wrong' way, we can remove the actual corner without meaning to.
        if jj==1 && closestWallDistsSorted(jj)==closestWallDistsSorted(jj+1)
            segRemoveInd( wallIDSort==closestWallOrder(jj+1) ) = false;
        end
        % Now actually remove the segments.
        ThSort       = ThSort( ~segRemoveInd );
        RSort        = RSort( ~segRemoveInd );
        wallIDSort   = wallIDSort( ~segRemoveInd );
    end
    % After having done all this, remove any remaining duplicate points at corners %
    %%% TW: I commented the two below lines out on 26/09/2023. I'm trying to make this run faster,
    %%%     so it becomes feasible to simulate BVCs at camera pixel resolution. The call to UNIQUE 
    %%%     was a big time sink with huge wall arrays, and I don't see why it is necessary - duplicate
    %%%     points should have been removed already in section 1C.
%     RTHUni = unique( [RSort, ThSort], 'rows', 'stable' );
%     RSort  = RTHUni(:,1);    ThSort = RTHUni(:,2);

    % 3d. Evaluate the BVC function at the r and phi coords of the wall/barrier, to get the response to the BVC %
    %     HERE IS THE ACTUAL BVC FUNCTION     %
    distFunc    = exp(    -((RSort-dist).^2)  ./  (2*(sigRad^2))    )   ./    sqrt(2*pi*(sigRad^2));   % This defines firing depending on distance from wall.
    bvcRespAng  = distFunc .* von_mises_PDF(ThSort, phi, kappa);   % This is the bvc response evaluated at angles correpsonding to each wall segment (segments defined by bins in map).
    
    % 3e. Before integrating to get the final, overall response, we need to weight the segment responses by the size of the angle subtended at the rat (i.e. close segments count more 
    % than far segments). Get the weights by sorting angles and doing DIFF, however, this measures the steps *between* wall coords, whereas function is evaluated *at* wall coords. 
    % So, to avoid circular asymmetry problems, do the DIFF in both directions, and take the mean.
    angDiffCCW             = [diff( ThSort ); circ_dist(ThSort(1),ThSort(end))];
    angDiffCW              = [circ_dist(ThSort(1),ThSort(end)); abs( flipud( diff( flipud(ThSort) ) ) )];
    % Before getting the means of these two, we also need to look for places where the visible wall segment switches from a far wall to a close 
    % barrier (or vice versa), when this happens the angle of the transition should be ignored (other one under-weights the close wall).
    RJumpInd               = (   abs( [diff(RSort);  RSort(1)-RSort(end)] ) ./ RSort   ) >= 0.2;
    angDiffCCW( RJumpInd )                          = nan;
    angDiffCW( [RJumpInd(end); RJumpInd(1:end-1)] ) = nan;
    % Now get the final weights by meaning the CW and CCW angles.
    %%% TW, 26/09/2023 - I changed 'circ_mean' to 'mean' in the line below, to save time, for
    %%% simulating BVCs at cam pix resolution. Can't think of a reason why it should really be
    %%% a circular mean, the angles are always tiny, none ever approach circle wrap-type size.
    angWeights             = mean( [angDiffCCW, angDiffCW], 2 );
      
    % 3f. Finally, the BVC response for this bin is the integrated sum of responses for each angle %
    bvcRespAngSortW              = bvcRespAng .* (angWeights./sum(angWeights));   % 'bvcRespAng' is already sorted circularly by angle.
    bvcMap( visY(ii), visX(ii) ) = nansum(bvcRespAngSortW);
    
    % debug plot %
    if 0 && visY(ii)==3 && visX(ii)==14
       [xPlot, yPlot] = pol2cart( ThSort, RSort );
       figure;  plot3( xPlot, yPlot, bvcRespAngSortW, 'b.-' );
       hold on;
       plot( 0, 0, 'ks' );
       set(gca, 'xgrid', 'on', 'ygrid', 'on' );
%        figure;  plot( 1:length(bvcRespAngSortW),  rad2deg(angWeights), 'kx' );
%        hold on
%        plot3( xPlot, yPlot, bvcRespAngSortW );
%        title(  bvcMap( visY(ii), visX(ii) ) );
%        figure;
%        plot( rad2deg(barrThSort), von_mises_PDF(barrThSort, phi, kappa) , 'kx');
    end

%     profile viewer
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (For 'model map' BVCs, function is finished here)
if ~simMode
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) For BVCs simulated with real position data, we now have convert 'bvcMap' (which
%     is essentially a LUT for BVC function value at visited locations) into a 
%     simulated spike train, then turn this into a rate map.
% 4a. Convert BVC function values into random poisson spike counts, scaled to matched 
%     requested peak rate.
bvcMap       = bvcMap( 2:end-1, 2:end-1 ); % Remove 1 pix-wide boundary where box edges were
bvcValsByPos = bvcMap(  sub2ind(mapSize-2,actY,actX)  );
bvcValsByPos = bvcValsByPos .* (simMeanRate/mean(bvcValsByPos,'omitnan'));  % Scale to desired mean rate
spksByPos    = poissrnd( bvcValsByPos./50 );                      % Convert to spikes per sample - note 50Hz samp rate hard coded here.       
% 4b. Standard binning and and smoothing of position and spikes to make rate map
%     NOTE 10 pixel bins and 5x5 boxcar smoothing hard coded here.
posBinned        = ceil([actY, actX] ./ 10);
posMap           = accumarray( posBinned, ones(size(actX))./50, [25 25] );
spkMap           = accumarray( posBinned, spksByPos, [25 25] );
visMask          = posMap~=0;
filtVis          = conv2(double(visMask),ones(5),'same');
filtPos          = conv2(posMap,ones(5),'same') ./ filtVis;
filtSpk          = conv2(spkMap,ones(5),'same') ./ filtVis;
bvcMap           = filtSpk./filtPos;
bvcMap(~visMask) = nan;













