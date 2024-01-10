function [Res, LUT_out] = generateAllBVCResponses(ModMaps,DataIn,varargin)
% Generate and save the complete BVC response map for every dataset.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Parameters %
prms.usePCModMaps         = 0;
prms.keepFullBVCRespArray = 0; % If 1, keeps r for every BVC tested, if 0 keeps only max r, then indices for best fit d, phi, SZ.
% BVC fit parameters
prms.BVCrBin         = 10;    % Pixels
prms.BVCrSmooth      = 5;  % Boxcar width in bins.
prms.BVCrAdaptSmooth = [];  % Adaptive smoothing kernel - if this is [], then defaults to using boxcat smooth.
prms.maxD4BarrFit    = 13;
% Vector map parameters.
prms.vMapSpatBin = 5;   % Sets spatial oversampling of vector map - should be an even divisor of prms.placeBin, so 1, 2, 5 or 10 (1=cam pix bins, 10=no oversampling relative to regular map).
prms.LUT         = [];  % Space-to-vector LUT, *for one dataset only*. Only used when this function called by 'generateShuffledBVCResponses'.
% General parameters - filtering + scaling data, which trials to analyse.
prms.minOccForEdge   = 50; % for path scaling
prms.speedFilter     = [2.5 100];
prms.nBslTrial       = 2;
prms.probeList       = {'barrier N-S', 'barrier W-E', 'CCE'}; % 
prms.boxExtent       = [250 280 320];  % Sets to sizes (in cam pix) to which data is scaled for [London_sq, Moser_sq, Circle] data. 
prms.useCircleTrans  = 0;  % This transforms all square maps into circles - needed to control whether 90-deg clustering of max phi is artifact of shape.
% Sub-sampling data to match path length, speed, etc.
prms.subSampMode   = 'MR'; % 'PL';
prms.subSampThrPL  = [1 1 1 1 1].*65;
prms.subSampThrSpd = [1 1 1 1 1].*7.57;
prms.subSampThrMR  = 2.64 ./ [2.68; 2.76; 2.64; 3.69];
prms.ageBins       = [16 18; 19 21; 22 25; 40 40]; % Some sub-samplings are age-bin specific, so require bins defined already.
% Stuff on how this function itself works.
prms.showTimer       = 1;
prms.makeCircRespMap = 0;
% Get some parameters from the model maps, to store for later reference %
% if prms.usePCModMaps
%     [prms.dList,prms.dList,prms.szList] = deal(ones(2,1)); 
% else
%     prms.dList      = ModMaps.d_list;
%     prms.phiList    = ModMaps.phi_list;
%     prms.szList     = ModMaps.sz_list;
% end
% prms.dList      = ModMaps.d_list;
% prms.phiList    = ModMaps.phi_list;
% prms.szList     = ModMaps.sz_list;
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
% Generate the 'envList' last, in case calling function has specified nBsl or probeList.
prms.envList    = cat(2, repmat( {'hp'}, 1, prms.nBslTrial ), prms.probeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup %%%
% Get the SCAn data and put it in a cell array (for PARFOR compatiblilty) % TW EDIT
if isempty( DataIn )  % If no data supplied by caller, ask scan GUI for selected data.
    SD=gss;
    for ii=1:length(SD.selData)
        bigDataArray{ii,1} = evalin('base',SD.selData{ii});
        bigDataArray{ii,2} = SD.selData{ii}; 
    end
else
    bigDataArray = DataIn; % (:,1);  % In this usage, pass the scan dataset and its 
end
% Convert the d and phi actual bins into an index (for rebinning corr outputs) %
if ~prms.usePCModMaps
    [~,~,dBinInd]   = unique( ModMaps.dKey );
    [~,~,phiBinInd] = unique( ModMaps.phiKey );
    [~,~,szBinInd]  = unique( ModMaps.szKey );
end
% Helpful to have an explicit flag for whether path data will be square or circle: this will
% be inverted if data undergoes shape transform. Matches format of prms.envList
isEnvSqFlag     = ~strcmp(prms.envList, 'CCE');
% Set up the (basic) table %
varList =   {
             'cellID',      'string'; ...
             'trialInd',    nan(1,length(prms.envList));
             'BVCResp',     cell(1, length(prms.envList)); ...
             'vectMap',     cell(1, length(prms.envList)); ...
                 };
if ~prms.keepFullBVCRespArray
varList2 =   {
             'SZMaxInd',     nan(1,length(prms.envList)); ...
             'phiMaxInd',    nan(1,length(prms.envList)); ...
             'BVCRespMax',   nan(1,length(prms.envList)); ...
             'dMaxInd',      nan(1,length(prms.envList)); ...
             };
varList = cat( 1, varList, varList2);
end
varList = varList';
Res = cell2table( varList(2,:) );
Res.Properties.VariableNames = varList(1,:);
Res.Properties.UserData = prms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~prms.usePCModMaps
    % For barrier trials, we only correlate with BVC models that are 1) orthogonal to barrier +-45deg,
    % 2) d-tuning <= half box width. Set up the indices for selecting relevant BVCs for this here, before main loop.
    % First get indices for which *values* of phi fit for each barr ori.
    % ***  {1}=Barr EW <BVCs NS>, {2} = Barr NS <BVCs EW>.  ***
    phiListInd{1} = abs(ModMaps.phi_list-deg2rad(90))<=deg2rad(45)  | abs(ModMaps.phi_list-deg2rad(270))<=deg2rad(45);                  % (ModMaps.phi_list >= 225*pi/180 & ModMaps.phi_list <= 315*pi/180)  | (ModMaps.phi_list >= 45*pi/180 & ModMaps.phi_list <= 135*pi/180 ); % could move out of loop
    phiListInd{2} = abs(ModMaps.phi_list-deg2rad(180))<=deg2rad(45) | ModMaps.phi_list<=deg2rad(45) |  ModMaps.phi_list>=deg2rad(315);  % (ModMaps.phi_list >= 135*pi/180 & ModMaps.phi_list <= 225*pi/180)  | ModMaps.phi_list >= 315*pi/180 | ModMaps.phi_list <= 45*pi/180;
    % Do the same for d as for phi ..
    dListInd      = ModMaps.d_list <= prms.maxD4BarrFit;
    % Convert 'list' index (i.e. which values are good) to 'key' index (which model maps have those good values).
    % 'Key' index values represent an index into the 'list' of values in phi_list and d_list. The 'key'
    % indices themselves are in index into the 3rd-dim stacking of all the model maps.
    for itOr=1:2
        phiKeyInd            = ismember( ModMaps.phiKey, find(phiListInd{itOr}) );
        dKeyInd              = ismember( ModMaps.dKey, find(dListInd) );
        % .. and merge d and phi to get final indices - these give which models to use for the two barrier oris.
        barrDPhiKeyInd{itOr} = dKeyInd & phiKeyInd;
    end
end

warning('OFF', 'MATLAB:table:RowsAddedExistingVars');
if prms.showTimer;  hWait = waitbar(0,'');   end
for ii=1:size(bigDataArray,1)
    if prms.showTimer;   waitbar(ii/size(bigDataArray,1), hWait, bigDataArray{ii,2});   end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1). Pre-process the data, find the relevant trials %
    % 1a. Scale the data.
    ScanData     = bigDataArray{ii,1};
    isMoserData=0;   if strcmp( ScanData.trials(1).user.environment, 'hp70' );  isMoserData=1;   end
	[ScanData,~] = devSub_findEdgesAndScalePath(ScanData,prms.minOccForEdge,prms.boxExtent(1+double(isMoserData)));                    % This function will internally test env type to make sure only squares are scaled.
    ScanData     = devSub_findEdgesAndScalePath_cylinder(ScanData, 10, 0.5, prms.boxExtent(3), 'radiusEstimateMethod', 'largestWellSampledRadius');  % This function will internally test env type to make sure only CCE are scaled.
    envsRunInDS  = cell(1,length(ScanData.trials));
    for jj=1:length(ScanData.trials)
        ScanData.trials(jj).window_x = 512;
        ScanData.trials(jj).window_y = 512;
        envsRunInDS{jj}              = ScanData.trials(jj).user.environment;
    end
    % 1b. Get the relevant trials: output from this section is list of trial indices that corresponding
    % to {'hp', 'hp' .. etc , 'probe1', 'probe2'}, where the number of pre-probe hps to use is defined in
    % 'prms.nBslTrial', and the followilengthng list of probes in 'prms.probeList'. If a particular trial type 
    % isn't in the dataset, trialsToUse(n) will be NaN.   
    trialsToUse = getStandardTrialSequence( envsRunInDS, prms );    % NOTE, 2017-12-14. Selection of standard trial sequence moved to a separate function, for replicable generation across main functions.
    
    % 1c. If requested, transform the data (i.e. transform squares into circles and vice versa).
    if prms.useCircleTrans
        ScanData    = devSub_transSq2Ci(ScanData, isMoserData);
        isEnvSqFlag = ~isEnvSqFlag;
    end

    % 1d. If requested, subsample the data to match speed, pathlength etc.
    if ~isempty( prms.subSampMode )
        ScanData = subSampleData( ScanData, prms );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% (2). Make rate maps.
    % 2a BVC fit rate maps%
    mapsBVC        = rates_main(ScanData,'bin',prms.BVCrBin,'adaptive_smooth', prms.BVCrAdaptSmooth,'smooth',prms.BVCrSmooth,'mode','rate','filt_speed',prms.speedFilter, 'trial_index', trialsToUse(~isnan(trialsToUse)) );    
    % 2b Vector map maps - unsmoothed, if requested, with smaller bin than other analyses.
%    mapsSpkUnSm    = rates_main(ScanData,'bin',prms.vMapSpatBin,'smooth',1,'mode','spike','filt_speed',prms.speedFilter);
%    mapsPosUnSm    = rates_main(ScanData,'bin',prms.vMapSpatBin,'smooth',1,'mode','pos','filt_speed',prms.speedFilter);
%    vMapOvSmpF     = prms.BVCrBin / prms.vMapSpatBin; % 'Spatial over-sample factor' for vMaps, needed to adjust padding and barr coords accordingly.

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% (3-4) By-trial loop in which we calculate BVC fits (3) and vector maps (4).
    % Pre-allocate cell array for all BVC response sets and vector maps in this scan dataset.
    [bvcFitForDS, vMapForDS] = deal(cell( length(ScanData.trials(1).cells), length(trialsToUse) ));
    LUT_out                  = cell(1,  length(trialsToUse) );  % Space-to-vector look-up-tables for each trial - need to pass back to calling function, when shuffling
    for itTr=1:length(trialsToUse)
        if isnan(trialsToUse(itTr));   continue;   end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before either analysis, need get the environment size and type. Need to take account of:
        % - Moser square boxes larger than london boxes
        % - some London rats odd barriers
        % - is the data transformed square->circle and vice versa?
        % i) Size of env (max map dimension).
        if ~isEnvSqFlag(itTr)
            mapSize = prms.boxExtent(3)/prms.BVCrBin;   % All circles (either 'real' or transformed squares) are one standard size by design.
        elseif isMoserData
            mapSize = prms.boxExtent(2)/prms.BVCrBin;   % All Moser squares are 280 pix side length
        else
            mapSize = prms.boxExtent(1)/prms.BVCrBin;   % All London squares are 250 pix side length
        end
        % ii) Get the right set of BVC models depending on experiment + use of shape transform
        if ~isEnvSqFlag(itTr)
            % i. First, if shape flag says circle, use circle
            EnvType = 'CCE'; 
        elseif strcmp( ScanData.trials( trialsToUse(itTr) ).user.environment, 'CCE' )
            % ii. If env label in raw data says 'CCE', but shape flag is square, data has been transformed - use baseline square.
            EnvType = 'hp';
        elseif strncmp( ScanData.trials( trialsToUse(itTr) ).user.environment, 'barrier', 7 )   &&   any( strcmp(ScanData.user.rat, {'123', '124', '137', '138'}) )
            % iii. Some London rats were run with shorter barriers, but these are not labeled specifcally in raw data.
            EnvType = [ ScanData.trials( trialsToUse(itTr) ).user.environment, ' short' ];
        else
            % iv. In all other cases, we can trust raw data label to identify correct model map.
            EnvType = ScanData.trials( trialsToUse(itTr) ).user.environment;
        end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 3. BVC fit analysis - correlate rate maps with model maps.
        if ~any(strcmp( {ModMaps.mapDefs(1:end).env}, EnvType )); continue; end % if we are using the PC model maps there are no model maps for the barrier trials so we need to skip these!
        BVCModMaps   = ModMaps.modelMaps{ strcmp( {ModMaps.mapDefs(1:end).env}, EnvType ) };
        % need to pad model maps in case they are place cells - maybe should
        % generate new version and save on disk)?
        if prms.usePCModMaps
            szModMaps = size(BVCModMaps) + [2 2 0 0];
            BVCModMaps = num2cell(BVCModMaps,[1 2]);
            BVCModMaps = cellfun( @(x) padarray(x, [1 1], NaN), BVCModMaps, 'UniformOutput', 0);
            BVCModMaps = reshape(cell2mat(BVCModMaps),szModMaps);
        end
        % 3a. Pad the ratemaps to fit barrier pixels at the edge (the edge of the array is the edge of the environment).
        mapsBVCForTr = cellfun( @(x) x(1:mapSize,1:mapSize), mapsBVC( trialsToUse(itTr),:), 'UniformOutput', 0 );  % Just get the TL corner - this is where the data should be after scaling. Different map sizes for squares and circles.
        mapsBVCForTr = cellfun( @(x) padarray(x, [1 1], NaN), mapsBVCForTr, 'UniformOutput', 0);
        % 3b. If trial is barrier, only correlate with subset of BVC models (orthogonal +-45, d<=half box width).
        % Indices for these set up before loop, see above.
%         if contains(EnvType,'N-S') || contains(EnvType,'V') 
%             BVCModMaps(:,:,~barrDPhiKeyInd{2}) = NaN;
%         elseif contains(EnvType,'W-E') || contains(EnvType,'H') 
%             BVCModMaps(:,:,~barrDPhiKeyInd{1}) = NaN;  
%         end
        % 3c. Now loop through cells, correlating each one with the whole model map set.
        for itCl=1:size( mapsBVC, 2)
            rateMap = mapsBVCForTr{1,itCl};
            % Get the pearson's r by hand so as to vectorise, this is the equation:
            % a = a - mean2(a);     b = b - mean2(b);
            % r = sum(sum(  a.*b  ))        /       (  sqrt( sum(sum(  a.*a  ))  *   sqrt( sum(sum(   b.*b   )) )    );
            A  = rateMap - nanmean(rateMap(:));                                      % Subtract samples sample mean. Dims of A here are 1,2
            B  = bsxfun(@minus, BVCModMaps, nanmean( nanmean(BVCModMaps,1), 2 ) );   % Dims of B are 1,2,3, mean for subtraction is calculated along third dim.
            AB = bsxfun(@times, A, B);
            sumAB =  nansum(nansum(AB,1),2);                    % This is a 1x1xnMod vector, one sum(AB) for each model.
            sqrAA = sqrt( nansum(nansum(A.^2,1),2) );           % This is a 1x1 integer
            sqrBB = sqrt( nansum(nansum(B.^2,1),2) );           % This is a 1x1xnMod vector, one sqrt(B.^2) for each model.
            R     = sumAB ./ bsxfun(@times, sqrAA, sqrBB);      % This is a 1x1xnMod vector, one pearsons-r for each model.
            if prms.usePCModMaps
                RArr(1,:,:) = squeeze(R); 
            elseif prms.keepFullBVCRespArray
                RArr  = accumarray( [dBinInd, phiBinInd, szBinInd], squeeze(R), [max(dBinInd), max(phiBinInd), max(szBinInd)], @nansum, nan );   % Reshape the r-vals into a ( 1:nPhiBin, 1:ndBin, 1:szBin ) array. (note that the 'nansum' isn't actually summing anything, there's just one value per bin, it's just a way of assigning values robust to NaNs.
            else
                [rMax, maxRInd] = max( R, [], 'omitnan' );
                dMax            = dBinInd( maxRInd );
                phiMax          = phiBinInd( maxRInd );
                szMax           = szBinInd( maxRInd );
                RArr            = [rMax, dMax, phiMax, szMax];
            end
            bvcFitForDS{ itCl, itTr } = RArr; 
            clear RArr
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 4. Vector map analysis.
        % 4a. Get the the unsmoothed maps in the correct spatial format: just the scaled data in TL corner ..
        % (1:mapsize, 1:mapsize), and then pad this with one line of NaN at each edge (for wall position).
%         cropAndPadMapFunc            = @(x) padarray( x( 1:(mapSize*vMapOvSmpF), 1:(mapSize*vMapOvSmpF) ), [1 1], NaN);
%         padSpkMaps                   = cellfun( cropAndPadMapFunc, mapsSpkUnSm(trialsToUse(itTr),:), 'uni', 0);
%         padPosMap                    = feval( cropAndPadMapFunc, mapsPosUnSm{trialsToUse(itTr)} );
%         % Check if a LUT has been supplied - if so get the one for the current trial.
%         if isempty(prms.LUT);   LUTForTr = [];   else;    LUTForTr = prms.LUT{itTr};   end
%         % Do a similar conversion for wall+barrier coords - here we need to split the code for squares (+straight barriers) and circles.
%         if isEnvSqFlag(itTr)
%             % 4b. Squares + straight barriers - make sure that wall remains defined as outer-most bin in map (i.e. 1 & vect_mapsize+2),
%             % barrier position multiplied by 'ovSmpF', but also adjusted for different wall bin size.
%             barrStr                        = ModMaps.mapDefs(   strcmp( {ModMaps.mapDefs(1:end).env}, EnvType )  ).barrStr;
%             RBWallInd                      = barrStr(:,5)==0 & barrStr>1;  % Right & bottom wall index. (Coords for T & L wall remain at 1)
%             barrStr(RBWallInd)             = (mapSize*vMapOvSmpF)+2;
%             barrInd                        = barrStr(:,5)>0;                        % barrInd = index for internal barriers.
%             barrStr(barrInd)               = ((barrStr(barrInd)-1) .* vMapOvSmpF) + 1;  % Barriers: -1 to remove existing corr for wall, *ovSmpF, then add 1 for wal again.
%            [vMapForDS(:,itTr),~,LUT_out{itTr}] = bvcTrVectMap_forLM( padPosMap, padSpkMaps, barrStr, [], 'nBinsDist', ceil(mapSize/2), 'vMapBinDist', vMapOvSmpF, 'vMapBinAng', 360/length(ModMaps.phi_list), 'LUT', LUTForTr );
%         else
%             % 4c. Circles - adjusted similar to straight internal barriers, but with an additional fixed offset  
%             % from TL of +0.5,as they are *centred* with respect to rows/cols of 'edge' bins. 
%             barrCirc                       = ModMaps.mapDefs(   strcmp( {ModMaps.mapDefs(1:end).env}, EnvType )  ).barrCirc;
%             barrCirc(1:2)                  = ((barrCirc(1:2)-1.5) .* vMapOvSmpF) + 1.5;
%             barrCirc(3)                    = ((barrCirc(3)-0.5) * vMapOvSmpF) + 0.5;
%            [vMapForDS(:,itTr),~,LUT_out{itTr}] = bvcTrVectMap_forLM( padPosMap, padSpkMaps, [], barrCirc, 'nBinsDist', ceil(mapSize/2), 'vMapBinDist', vMapOvSmpF, 'vMapBinAng', 360/length(ModMaps.phi_list), 'LUT', LUTForTr );
%         end

        
        % END for by-trial loop.
    end        


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) Assign results for this dataset to output table % 
    cellID = cell(length(ScanData.trials(1).cells), 1);
    for jj=1:length(cellID);   cellID{jj} = [bigDataArray{ii,2}, ' t', num2str(ScanData.trials(1).cells(jj).tet), 'c', num2str(ScanData.trials(1).cells(jj).cellnum)];   end
    rowInd = 1:length(cellID);
    if ii>1;   rowInd = rowInd + (size(Res,1));   end      % Need to take account that Res is initialised at one row long already.
    Res.cellID( rowInd, 1 )                       = cellID;
    Res.trialInd( rowInd, 1:length(trialsToUse) ) = repmat( trialsToUse, length(rowInd), 1 );
    if prms.keepFullBVCRespArray
        Res.BVCResp( rowInd, 1:length(trialsToUse) )  = bvcFitForDS;
    else
        bvcFitForDS( cellfun(@isempty, bvcFitForDS ) )   = {[nan nan nan nan]};
        Res.BVCRespMax( rowInd, 1:length(trialsToUse) )  = cellfun( @(x,y) x(1), bvcFitForDS );
        Res.dMaxInd( rowInd, 1:length(trialsToUse) )     = cellfun( @(x,y) x(2), bvcFitForDS );
        Res.phiMaxInd( rowInd, 1:length(trialsToUse) )   = cellfun( @(x,y) x(3), bvcFitForDS );
        Res.SZMaxInd( rowInd, 1:length(trialsToUse) )    = cellfun( @(x,y) x(4), bvcFitForDS );
    end
    Res.vectMap( rowInd, 1:length(trialsToUse) )  = vMapForDS;
        
    
end
if prms.showTimer;   close(hWait);   end
warning('ON', 'MATLAB:table:RowsAddedExistingVars');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (6) If requested, re-shape the BVC response map to a circular display %
% % Get the coorindate grids for the transform %
% if prms.makeCircRespMap
%     [phiBinsArr,dBinsArr]    = meshgrid( prms.phiList, prms.dList );
%     [polBinsInX, polBinsInY] = pol2cart( phiBinsArr, dBinsArr+1 ); 
%     [resampX, resampY]       = meshgrid( -max(prms.dList):max(prms.dList) );
%     % And run it %
%     Res.BVCRespCirc = cell( size(Res.BVCResp) );
%     for ii=1:numel(Res.BVCResp)
%         RArr                = flipud( Res.BVCResp{ii} );
%         if isempty(RArr);   continue;   end
%         resampR             = griddata( polBinsInX(:), polBinsInY(:), RArr(:), resampX(:), resampY(:) );
%         Res.BVCRespCirc{ii} = reshape( resampR, size(resampX) );
%     end
% 
% end
%     
    
    
    













