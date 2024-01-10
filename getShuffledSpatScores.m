function [Res] = getShuffledSpatScores(DataIn, varargin)
% Spike shuffled thresholds for directionality for single rate maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shuffling parameters %
prms.fastMode       = 1;  % If fastMode=1, skips measures based on adaptive smoothed maps.
prms.shufNStepsSpat = 1000;
prms.shufMinOffset  = 20;
% Rate map parameters %
prms.placeBin       = 10;
prms.placeSmooth    = 5;
prms.adaptSmooth    = 200;    % 200 is what Langston 2010 used for SI calculation.
prms.dirBin         = 6;
prms.dirSmooth      = 5;
prms.speedFilter    = [2.5 100];
% Border score.
prms.bdScRateThr  = 0.2;  % Proportion of max rate to define field. 0.3 in Solstad 2008, 0.2 in Bjerknes 2014.
prms.bdScSizeThr  = 200;  % Field area threshold in cm2 - function requires bin size in sq cm. For both Moser papers = 200.
% Specification of 'standard' trial sequence.
prms.nBslTrial  = 2;
prms.probeList  = {'CCE'};  %  {'barrier N-S', 'barrier W-E' , 'CCE'}; %, 'net'};
prms.envList    = cat(2, repmat( {'hp'}, 1, prms.nBslTrial ), prms.probeList);
prms.useCircleTrans = 1;
% Scaling data to fit standard square %
prms.minOccForEdge = 50;
prms.boxExtent     = [250, 280];
% Exclusion and filtering of data %
prms.ratExcludeList           = {'139','158','218'};
prms.cellsMustFireInAllTrials = 1;
prms.nSpkDefForFiring         = 100;
% Sub-sampling data to match path length, speed, etc.
prms.subSampMode   = [];
prms.subSampThrPL  = [1 1 1 1 1].*65;
prms.subSampThrSpd = [1 1 1 1 1].*7.57;
prms.subSampThrMR  = [0.7103; 1.0000; 0.9307; 0.3623]; % format is (nAgeBin,1) - I haven't got these numbers by trial type yet.
prms.ageBins       = [16 18; 19 21; 22 25; 40 40];     % Some sub-samplings are age-bin specific, so require bins defined already.
% The most trials in any dataset. Needed for pre-allocation. %
prms.maxNTrial    = length(prms.envList);     
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
[bdScTemp, SITemp, SIAdSmTemp, datasetNamesTemp, bdScAdSmTemp, RVTemp] = deal( cell( 1,size(bigDataArray,1) ) );   % Pre-allocate temp cell arrays for results (cannot combine into structure until after PARFOR).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop for getting shuffled data scores %
%tic
parfor ii=1:length(bigDataArray)  %%% PAR
% for ii = 1
    data=bigDataArray{ii};
    if strcmp(data.trials(1).user.environment,'hp70');   isMoserData=true;   else   isMoserData=false;   end
    mapSize     = prms.boxExtent( 1+double(isMoserData) ) / prms.placeBin;    % #ok<PFBNS

    % Declare temporary variables for storing scores, for this dataset (PARFOR nonsense).   
    [bdScTempDS, SITempDS, SIAdSmTempDS,  bdScAdSmTempDS, RVTempDS ] = deal(nan( length(prms.envList), length(data.trials(1).cells), prms.shufNStepsSpat));   % Pre-allocate temp arrays for this dataset.

    % Scale the data (fit squares to be 250 x 250 pix exactly) %
    [data,SF] = devSub_findEdgesAndScalePath( data, prms.minOccForEdge, prms.boxExtent( 1+double(isMoserData) ) );
    % If requested, transform to square to circle and vice-versa.
    if prms.useCircleTrans
        data        = devSub_transSq2Ci(data, isMoserData);
    end
    % If it is moser data, it will need the ppm (and therefore speed) correcting %
    if isMoserData
        for jj=1:length(data.trials)
            trueOriginalPPM       = 400 / mean( SF(jj,:) ); % Take the mean across x and y. If you have had to double the path coords to acheive 400ppm, then your original ppm was actaully 200 (same distance less pixels). 
            data.trials(jj).speed = data.trials(jj).speed .* (300/trueOriginalPPM);  % The speed has been calcualted assuming ppm=300, it should be rescaled accordingly to fix this. If ppm was too large, speed will have been too small (as speed proportional to m/pix)
            data.trials(jj).ppm   = 400;   % Set this so that path length is calcualted correctly, below.
        end      
    end
    
    % If requested, sub-sample the data
    if ~isempty( prms.subSampMode )
       data = subSampleData( data, prms ); 
    end
    
    % Get the list of trials to analyse %
    envList     = cell(1,length(data.trials));
    for jj=1:length(data.trials);   envList{jj}=data.trials(jj).user.environment;    end
    trialsToUse = getStandardTrialSequence( envList, prms );

    
    for jj=1:length(trialsToUse)  %% jj = iterator for trial (but note that the first action here is to reduce the data to one trial only, so after that jj only used for output assignment).
        
        if isnan( trialsToUse(jj) );  continue;  end
        
        trData=data;
        trData.trials = trData.trials(  trialsToUse(jj)  );
        
        offsets = linspace( prms.shufMinOffset, trData.trials.dur-prms.shufMinOffset, prms.shufNStepsSpat);
        
        if strcmp( data.trials( trialsToUse(jj) ).user.environment, 'CCE' ) && ~prms.useCircleTrans
            isSqEnv = 0;
        elseif ~strcmp( data.trials( trialsToUse(jj) ).user.environment, 'CCE' ) && prms.useCircleTrans
            isSqEnv = 0;
        else
            isSqEnv = 1;
        end

        % Also need to loop through *Each indivdual offset* %
        for kk=1:length(offsets) %% kk = iterator for offset
            
            dataShuf = spk_randomise(trData,'fixedWrap',offsets(kk));
            
            % Make maps (just for 1 trial, one offset %
            maps                     = rates_main(dataShuf,'bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilter  );
            mapsPos                  = rates_main(dataShuf,'bin',prms.placeBin,'smooth',prms.placeSmooth, 'mode', 'pos','filt_speed',prms.speedFilter );
            mapsDir                  = rates_main(dataShuf,'bin_dir',prms.dirBin,'smooth',prms.dirSmooth,'space','dir','filt_speed',prms.speedFilter );
            if ~prms.fastMode
                [mapsAdSm, mapsPosAdSm]  = rates_main(dataShuf,'bin',prms.placeBin,'adaptive_smooth',prms.adaptSmooth,'filt_speed',prms.speedFilter  );
            end

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get the scores for each cell
            for mm=1:length(trData.trials.cells)
                % Boxcar maps %
                map                  = maps{1,mm};
                mapPos               = mapsPos{1};
                mapDir               = mapsDir{1,mm};
                if isSqEnv
                    bdScTempDS(jj,mm,kk) = calculate_BorderScore( map( 1:mapSize, 1:mapSize ), ((prms.placeBin/400)*100)^2, 'rateThr', prms.bdScRateThr, 'sizeThr', prms.bdScSizeThr );
                end
                SITempDS(jj,mm,kk)   = map_skaggsinfo( map, mapPos  );
                RVTempDS(jj,mm,kk)   = dir_rayleighvector( mapDir ); 


                % Adaptive smoothed maps %
                if prms.fastMode
                    bdScAdSmTempDS(jj,mm,kk)    = nan;
                    SIAdSmTempDS(jj,mm,kk)      = nan;
                else
                    map                      = mapsAdSm{1,mm};
                    mapPos                   = mapsPosAdSm{1};
                    if isSqEnv
                        bdScAdSmTempDS(jj,mm,kk) = calculate_BorderScore( map( 1:mapSize, 1:mapSize ), ((prms.placeBin/400)*100)^2, 'rateThr', prms.bdScRateThr, 'sizeThr', prms.bdScSizeThr );
                    end
                    SIAdSmTempDS(jj,mm,kk)   = map_skaggsinfo( map, mapPos  );
                end

            end

        end
        
        % Transfer the temp arrays for this dataset to the overall cell array %
        bdScTemp{ii}           = bdScTempDS;
        SITemp{ii}             = SITempDS;
        bdScAdSmTemp{ii}       = bdScAdSmTempDS;
        SIAdSmTemp{ii}         = SIAdSmTempDS;
        RVTemp{ii}             = RVTempDS;

        datasetNamesTemp{ii}   = bigNameArray{ii};

    end
end
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-organise the output as a table, for ease of comparison with main data %
% Set up results table %
scoreDum = cell(1,prms.maxNTrial);
for ii=1 % For neatness, the pre-allocation variable list is hidden in this folded (dummy) FOR loop. 
    varList =   {
         
                 'cellID',           'string'; ...
                 'age',              nan; ...
                 'ShufborderScore',  scoreDum; ...
                 'ShufSI',           scoreDum; ...
                 'ShufRV',           scoreDum; ...

                 'ShufborderScoreAdSm',  scoreDum; ...
                 'ShufSIAdSm',           scoreDum; ...


                     };
    varList = varList';
end
dummyRow = cell2table( varList(2,:) );
Res = dummyRow;
Res.Properties.VariableNames = varList(1,:);
Res.Properties.UserData = prms;

% Loop though results arrays from above and assign data to rows in table %
cellCount=1;
for ii=1:length(datasetNamesTemp)  %% ii = iterator for dataset
    
    for jj=1:size(bdScTemp{ii},2)  %% jj = iterator for cell
        
        Res(cellCount,:) = dummyRow;
        
        % Assign metadata for Cell %
        data=bigDataArray{ii};
        Res.cellID(cellCount,1) = {[datasetNamesTemp{ii} ' t' num2str(data.trials(1).cells(jj).tet) 'c' num2str(data.trials(1).cells(jj).cellnum)]};
        Res.age(cellCount,1)    = str2double( data.user.age );
        
        % Assign data by trial %
        for kk=1:size(bdScTemp{ii},1)  %% kk=iterator for trial
            
            Res.ShufborderScore(cellCount,kk) = { squeeze( bdScTemp{ii}(kk,jj,:)) };
            Res.ShufSI(cellCount,kk)          = { squeeze( SITemp{ii}(kk,jj,:)) };
            Res.ShufRV(cellCount,kk)          = { squeeze( RVTemp{ii}(kk,jj,:)) };

            Res.ShufborderScoreAdSm(cellCount,kk)  = { squeeze( bdScAdSmTemp{ii}(kk,jj,:)) };
            Res.ShufSIAdSm(cellCount,kk)           = { squeeze( SIAdSmTemp{ii}(kk,jj,:)) };
            
        end
        
        % Bump the cell count %
        cellCount = cellCount + 1;
    end
    
end
            








    
