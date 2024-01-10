function [Res] = generateShuffledBVCResponses(ModMaps, DataIn, varargin)
% Generate shuffled BVC fit and vector response maps, and save for each rate map:
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analysis parameters %%%
prms.usePCModMaps         = 0;
prms.keepFullBVCRespArray = 0; % If 1, keeps r for every BVC tested, if 0 keeps only max r, then indices for best fit d, phi, SZ.
% Shuffling parameters %
prms.shufMinOffset  = 20;           % In sec
% prms.shufNSteps = 200;           % 200
prms.ageBinsShuf    = [13 18; 19 21; 22 25; 26 36; 40 40];
%prms.shufNStepsBVC  = {ceil([124,147,293,0,105]./20),ceil([480,339,350,176,260]./20)}; % SUB data and Moser data for canonical agebins [16 18; 19 21; 22 25; 26 36; 40 40]
prms.shufNStepsBVC  = {   [1000,1000,1000,0,1000]./1,  [1000,1000,1000,1000,1000]./1  };
% Params for 'generateAllBVCResponses'. Note that data-sensitive binning, filtering otocptions are not
% here, these should only be manipulated within that function.
prms.nBslTrial      = 2;
prms.probeList      = {'CCE'};  %  'barrier N-S', 'barrier W-E', , 'net'};
prms.envList        = cat(2, repmat( {'hp'}, 1, prms.nBslTrial ), prms.probeList);
prms.showTimer      = 0;
prms.LUT            = [];  % LUT for vMaps initialised to [], then after first shuffle run per dataset, will get this returned, then fed to calls in following shuffles.
% Sub-sampling data to match path length, speed, etc.
prms.subSampMode   = 'MR';
prms.subSampThrPL  = [1 1 1 1 1].*65;
prms.subSampThrSpd = [1 1 1 1 1].*7.57;
prms.subSampThrMR  = [0.7103; 1.0000; 0.9307; 0.3623]; % format is (nAgeBin,1) - I haven't got these numbers by trial type yet.
% Get some parameters from the model maps, to store for later reference %
% if prms.usePCModMaps
%     [prms.dList,prms.dList] = deal(ones(2,1));
% else
%     prms.dList      = ModMaps.d_list;
%     prms.phiList    = ModMaps.phi_list;
% end
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% --------------------------------------------------------------    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the SCAn data and put it in a cell array (for PARFOR compatiblilty) % TW EDIT
if isempty( DataIn )  % If no data supplied by caller, ask scan GUI for selected data.
    SD=gss;
    for ii=1:length(SD.selData)
        bigDataArray{ii} = evalin('base',SD.selData{ii});
        bigNameArray{ii} = SD.selData{ii}; 
    end
else
    bigDataArray = DataIn(:,1);  % In this usage, pass the scan dataset and its 
	bigNameArray = DataIn(:,2);  % name as a {1,2} cell, to the prms struct field 'dataIn'
end
% Pre-allocate temp cell arrays to hold results for each dataset (cannot combine into structure until after PARFOR). %
[BVCRespByDS, vMapByDS, BVCRespMaxByDS,BVCRespMaxDByDS,BVCRespMaxPhiByDS,BVCRespMaxSZByDS,vMapPk2MnRByDS,vMapMaxDByDS,vMapMaxPhiByDS] = ...
    deal(  cell(1, size(bigDataArray,1) )  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop for getting shuffled data scores %
%WaitMessage = parfor_wait( size(bigDataArray,1) );
nCores = feature('numcores');
parfor itDS=1:length(bigDataArray)  %%% PAR START
% for itDS= 1:length(bigDataArray)
   % WaitMessage.Send;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set-up for shuffling %%%
    data=bigDataArray{itDS};
    
    % Number of shuffles per trial depends on brain region and age bin, such that all of these region/age groups come
    % out with roughly the same number of total shuffles. Here we detect the brain region and age, and set 'shufNSteps' appropriately.
    if strcmp(data.trials(1).user.environment,'hp70')
        brainRegionInd = 2;
        nShufInd = str2double(data.user.age) >= prms.ageBinsShuf(:,1) & str2double(data.user.age) <= prms.ageBinsShuf(:,2);
    else
        brainRegionInd = 1;
        if str2double(data.user.age) == 32
            nShufInd = 40 >= prms.ageBinsShuf(:,1) & 40 <= prms.ageBinsShuf(:,2); % need to hard code canonical adult age here for SUB data
        else
            nShufInd = str2double(data.user.age) >= prms.ageBinsShuf(:,1) & str2double(data.user.age) <= prms.ageBinsShuf(:,2);
        end
    end
    shufNSteps = prms.shufNStepsBVC{brainRegionInd}(nShufInd); % this selects dataset specific n of shuffles (depending on brain region and age)
    
    % Get the offsets (i.e. absolute time shifts in sec): these have to be organised by trial, as the trials can be different lengths %
    offsets = cell(1,length(data.trials));
    for jj=1:length(data.trials)
        offsets{jj} = linspace( prms.shufMinOffset, data.trials(jj).dur-prms.shufMinOffset, shufNSteps );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Now loop through each indivdual offset %%%
    % Pre-allocate temp cell array for this dataset: each cell in the array will hold the data for one rate map, for iterations of shuffle.
    [BVCRespAllShufs, vMapAllShufs] = deal(cell( length(data.trials(1).cells), length(prms.envList) ));
    
    prms_cp = prms;    % The way we supply the vMap LUT for a dataset to 'generateAllBVCResponses' (i.e. in a prms struct) won't work with parfor (can't mod a broadcast var), so need to take a 'local loop' copy first. 
    for itSh=1:shufNSteps     %% itSh = iterator for offset

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the data for this (one) shuffle): time-shift the spike trains and run main analysis function (generateAllBVCResponses).
        offsetsForLoop    = cat(1, cellfun( @(x,y) x(y), offsets, repmat({itSh},size(offsets)), 'UniformOutput', 1));
        dataShuf          = spk_randomise(data,'fixedWrap',offsetsForLoop);
        if itSh==1
            [ShufDataTemp,prms_cp.LUT] = generateAllBVCResponses( ModMaps, {dataShuf, bigNameArray{itDS}}, prms_cp );
        else
            ShufDataTemp               = generateAllBVCResponses( ModMaps, {dataShuf, bigNameArray{itDS}}, prms_cp );
        end

        
        % Following this, take the BVC response maps (from 1 shuffle iteration) and add to compiled array of BVC responses from all shuffles
        % Note the concatenation dimension is 4, to take account of possible multiple sigma zero, which when exists is 3rd dim of 'BVCResp'
        BVCRespAllShufs  = cellfun(@(x,y) cat(4,x,y), BVCRespAllShufs, ShufDataTemp.BVCResp, 'UniformOutput', 0);
%         vMapAllShufs     = cellfun(@(x,y) cat(4,x,y), vMapAllShufs, ShufDataTemp.vectMap, 'UniformOutput', 0);
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Transfer the shuffled data entries (maps & derived scores) for this dataset to the PARFOR-level cell array %
    %%% This depends on whether ouput from 'generateAllBVCResponses' is
    %%% full R array, or just R Max, d, phi and SZ for max (although,
    %%% former is just legacy mode, really, when shuffling should save time
    %%% and memory by not keeping full array in the first place).
    if prms.keepFullBVCRespArray
        % This is legacy - not sure why this would be used now
        [BVCRespMaxByDS{itDS}, BVCRespMaxDByDS{itDS}, BVCRespMaxPhiByDS{itDS}, BVCRespMaxSZByDS{itDS}]  = cellfun( @nanmax4D, BVCRespAllShufs, 'UniformOutput', 0 );
    else
        % This is the normal case
        [BVCRespMaxByDS{itDS}, BVCRespMaxDByDS{itDS}, BVCRespMaxPhiByDS{itDS}, BVCRespMaxSZByDS{itDS}] = deal( cell( size( BVCRespAllShufs ) ) );
        for itRM=1:numel( BVCRespAllShufs )
            if isempty( BVCRespAllShufs{itRM} );   continue;   end
            BVCRespMaxByDS{itDS}{itRM}    = squeeze( BVCRespAllShufs{itRM}(:,1,:,:) );
            BVCRespMaxDByDS{itDS}{itRM}   = squeeze( BVCRespAllShufs{itRM}(:,2,:,:) );
            BVCRespMaxPhiByDS{itDS}{itRM} = squeeze( BVCRespAllShufs{itRM}(:,3,:,:) );
            BVCRespMaxSZByDS{itDS}{itRM}  = squeeze( BVCRespAllShufs{itRM}(:,4,:,:) );
        end
    end
    
%     vMapAllShufs                                        = cellfun( @(x) x./mean(x,[1 2],'omitnan'), vMapAllShufs, 'UniformOutput', 0 );  % Mean normalise vMaps first thing.
%     vMapPk2MnRByDS{itDS}                                = cellfun( @(x) squeeze( max(x,[],[1 2],'omitnan') ), vMapAllShufs, 'UniformOutput', 0 );
%    [~, vMapMaxDByDS{itDS}, vMapMaxPhiByDS{itDS}]        = cellfun( @nanmax2D, vMapAllShufs, 'UniformOutput', 0 );
    
    % Store also indivdual shuffles for BVC responses (for vector-matched secondary classification).
    % This comes last, as we delete barrier trials, and convert 4dp double to int16, for space saving.
    BVCRespAllShufs        = cellfun(@(x) int16(round(x*10000)), BVCRespAllShufs, 'Uniformoutput', 0);
%     vMapAllShufs           = cellfun(@(x) int16(round(x*10000)), vMapAllShufs, 'Uniformoutput', 0);
    barrTrInd              = contains( prms.envList, 'barrier');
    if ~isempty(barrTrInd)
        BVCRespAllShufs(:,barrTrInd) = {[]};
%         vMapAllShufs(:,barrTrInd)    = {[]};
    end
%     BVCRespByDS{itDS}      = BVCRespAllShufs;
%     vMapByDS{itDS}         = vMapAllShufs;
    % Store DS and trial metadata.
    TrialIndByDS{itDS}     = ShufDataTemp.trialInd( 1, : );
    DatasetNamesByDS{itDS} = bigNameArray{itDS};
 
end
%WaitMessage.Destroy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-organise the output as a table, for ease of comparison with main data %
% Set up results table %
nTrialCols = prms.nBslTrial + length(prms.probeList);
scoreDum = cell(1, nTrialCols);
for ii=1 % For neatness, the pre-allocation variable list is hidden in this folded (dummy) FOR loop. 
    varList =   {
        
                 'cellID',      'string'; ...
                 'dataset',     'string'; ...     
                 'trialInd',     nan(1,nTrialCols); ...
                 
                 'ShufBVCResp',         scoreDum; ...
                 'ShufVMap',            scoreDum; ...

                 'ShufBVCRespMax',      scoreDum; ...
                 'ShufBVCRespMaxD',     scoreDum; ...
                 'ShufBVCRespMaxPhi',   scoreDum; ...
                 'ShufBVCRespMaxSZ',    scoreDum; ...
                 
                 'ShufvectMap_pk2MnR',  scoreDum; ...
                 'ShufvectPhiMax',      scoreDum; ...
                 'ShufvectDMax',        scoreDum; ...

                     };
    varList = varList';
end
dummyRow = cell2table( varList(2,:) );
Res = dummyRow;
Res.Properties.VariableNames = varList(1,:);
Res.Properties.UserData = prms;

% Loop though results arrays from above and assign data to rows in table %
nTr       = nTrialCols;
cellCount = 1;
for ii=1:length(DatasetNamesByDS)  %% ii = iterator for dataset
    
    for jj=1:size(BVCRespMaxByDS{ii},1)  %% jj = iterator for cell
        
        Res(cellCount,:) = dummyRow;
        
        % Assign metadata for Cell %
        data=bigDataArray{ii};
        Res.cellID(cellCount,1)       = {[DatasetNamesByDS{ii} ' t' num2str(data.trials(1).cells(jj).tet) 'c' num2str(data.trials(1).cells(jj).cellnum)]};
        Res.dataset(cellCount,1)      = DatasetNamesByDS(ii);
        Res.trialInd(cellCount,1:nTr) = TrialIndByDS{ii};
        
%         Res.ShufBVCResp(cellCount,1:nTr)        = BVCRespByDS{ii}( jj, : );
%         Res.ShufVMap(cellCount,1:nTr)           = vMapByDS{ii}( jj, : );
        Res.ShufBVCRespMax(cellCount,1:nTr)     = BVCRespMaxByDS{ii}( jj, : );
        Res.ShufBVCRespMaxD(cellCount,1:nTr)    = BVCRespMaxDByDS{ii}( jj, : );
        Res.ShufBVCRespMaxPhi(cellCount,1:nTr)  = BVCRespMaxPhiByDS{ii}( jj, : );
        Res.ShufBVCRespMaxSZ(cellCount,1:nTr)   = BVCRespMaxSZByDS{ii}( jj, : );

%         Res.ShufvectMap_pk2MnR(cellCount,1:nTr)  = vMapPk2MnRByDS{ii}( jj, : );
%         Res.ShufvectPhiMax(cellCount,1:nTr)      = vMapMaxPhiByDS{ii}( jj, : );
%         Res.ShufvectDMax(cellCount,1:nTr)        = vMapMaxDByDS{ii}( jj, : );
        
        % Bump the cell count %
        cellCount = cellCount + 1;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxVal, rInd, cInd, zInd] = nanmax4D( in )
% For a 4-D input, get the max in each 4-th dim slice, i.e. max(a(:,:,:,n)), and
% also return the row, col, and z-dim indices for each of these max values.
[maxVal,lInd]          = max( in, [], [1 2 3], 'linear', 'omitnan' );
[rInd,cInd,zInd,d4dum] = ind2sub( size(in), squeeze(lInd) );  % d4dum seems to be necessary to get ind2sub to work across all 4 dims.
maxVal                 = squeeze(maxVal);


            






    
