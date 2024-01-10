function [Res,rThrAll,rThr2] = shufThrCellClassifier(Res,ResShuf, varargin)
% As it says in the function name .. generates 'isBVC' field in Res struct.
%
%   Res                    = shufThrCellClassifier(Res,ResShuf,varargin)
%   [Res,thrVals]          = shufThrCellClassifier(Res,ResShuf,varargin)
%
% ResShuf can be set to [] is not needed, for example when using a fixed threshold.

% ANALYSIS PARAMETERS:

% Thresholding (this function) %
prms.ageBins            =  [16 18; 19 21; 22 25; 40 40];  % [15 17; 18 19; 20 21; 22 25; 32 32];
prms.thrScore           = {'BVCRespMax_VFS','SI'};  %'BVCRespMax'; %   '
prms.thrTrials          = 1:2;   % 1 - 5, single value or vector, following templateBSL, nonTemplateBSL, Barrier N-S, Barrier W-E, CCE. 
prms.thrTrialsShuf      = 1:2;
prms.thrTrialAv         = 'max'; %  'min'; % 'mean';  % How to get an average test value if length(prms.thrTrials)>1. One can also think that 'min' and 'max' are equivalent to ALL and ANY (respectively), and in the case of dual threshold classification, it explicitly works like this: 'min'=all must pass, 'max'=any can pass.
prms.thrPrc             = [99 75]; % 99; %    % Percentile of shuffled data r to use as threshold
prms.thrByMap           = [1 0];    % If 1, then thresholds are calculated for individual maps, and BVCs defined on if their scores crosses their self-generated threshold 
prms.fixedThrs          = [];   % Can supply fixed thrs instead of using shuffled data, either one single thr or one per age bin. Can also be for multiple measures. Format = (nMeasure,1) or (nMeasure,nAgeBin).
prms.outputFieldName    = 'isBVC';   % Allows the caller to custom name the output logical classifier field - useful if you want to do 'isBorderCell', etc.
prms.excPCFits          = 0;
% Low rate/poor path filtering (this function) %
prms.meanRateThr        = 0.2;
prms.nBinVisThrAbs      = [625 784 784].*0.8;  % 627 specific for Moser data: 0.8 * (28^2)
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for itAB=1:2:length(varargin);   prms.(varargin{itAB}) = varargin{itAB+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for itAB=1:length(f);   prms.(f{itAB}) = s.(f{itAB});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
% Input format check - to avoid having to change calling code after I changed this function to allow multiple thresholds.
if ischar( prms.thrScore )
    prms.thrScore = {prms.thrScore};  
end
% If fixed thrs supplied for only one age or brain region, repmat out to (nMeasure,nAgeBin,2) for ease of indexing.
if size( prms.fixedThrs, 2 ) == 1
    prms.fixedThrs = repmat( prms.fixedThrs, 1, size(prms.ageBins,1), 1);
end
if size( prms.fixedThrs, 3 ) == 1
    prms.fixedThrs = repmat( prms.fixedThrs, 1, 1, 2);
end
% If dealing with transformed data, need to switch thresholds for bin sampling accordingly.
if Res.Properties.UserData.useCircleTrans
    prms.nBinVisThrAbs = prms.nBinVisThrAbs( [3 2 1] ); % Swap Sq and Ci, Moser stays the same.
end
% Merge shuffled table input with actual data table.
if ~isempty(ResShuf) 
    thrScoreNames = unique(prms.thrScore);
    for itSc=1:length(thrScoreNames)
        ShufScores{itSc}  = ['Shuf' thrScoreNames{itSc}];  
    end
    % Check that none of the shuffled scores are in the Res struct already (but why should they be?)
    fNames       = fieldnames(Res);
    shufFieldInd = cellfun(@(x) any(strcmp(x,fNames)), ShufScores);
    ShufScores   = ShufScores(~shufFieldInd);
    % If the Voronoi field shuffle is a specified thresholding score, make a copy of the normal
    % 'BVCRespMax' field, named 'BVCRespMax_VFS', so that the naming logic for shuffled and data 
    % scores is consistent (will remove this at the end).
    if any(strcmp(ShufScores,'ShufBVCRespMax_VFS'))
        Res.BVCRespMax_VFS = Res.BVCRespMax;
    end
    % Join the required shuffled scores to the main Res structure.
    if ~isempty(ShufScores)
        Res = join(Res, ResShuf, 'Keys', 'cellID', 'RightVariables', ShufScores);    % If you put Res as the first argumnent then it shouldn't change the row order.
    end
%     % envList in Res and ResShuf might be different (e.g. don't shuffle barrs to save time). Need to match indices. 
%     prms.thrTrialsShuf = find( strcmp( Res.Properties.UserData.envList(prms.thrTrialsShuf(1)), ResShuf.Properties.UserData.envList ) );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a filter index of when a cell is 'well-sampled': is the mean rate and position smapling sufficient? %
wellSampledInd             = Res.meanRate >= prms.meanRateThr;
nBinVisThr                 = ones(size(Res.nBinVis)) .* prms.nBinVisThrAbs(1);  % Default coverage threshold
nBinVisThr(Res.area==2, :) = prms.nBinVisThrAbs(2);                       % Moser coverage threshold
nBinVisThr(    :,       3) = prms.nBinVisThrAbs(3);                       % CCE coverage threshold
wellSampledInd             = wellSampledInd & Res.nBinVis>=nBinVisThr;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find thresholds. Two ways of doing this, (1) map-by-map (assumes you have done enough shuffles to 
%%% make the threshold meaningful, (2) pooled-cell thresholds, grouped by age and area. 
%%% Function allows mix and match of these two approaches (so e.g. BVC response is by map, SI is by population).
%%% To accomodate this, call each classification method in a sub function, for the appropriate classification measure.

thrByMapMsInd                          = logical( prms.thrByMap );
BVCIndByTrial                          = true( size(Res,1), length(prms.thrTrials) );
rThrAll                                = nan(length(prms.thrScore),size(prms.ageBins,1),2);

for itMs=1:length(prms.thrScore)
 
    if thrByMapMsInd( itMs )
        % By map measure(s)
        prmsCopy_byMap                            = prms;
        prmsCopy_byMap.thrScore                   = prms.thrScore( itMs );
        prmsCopy_byMap.thrPrc                     = prms.thrPrc( itMs );
        [BVCIndByTrial_MsTemp, thrTmpPop, thrTmpMap] = classifyByMapSubFunc( Res, wellSampledInd, prmsCopy_byMap );
        rThrAll( itMs, :, : )                     = thrTmpPop;
        Res.(['thr_' num2str(prms.thrPrc(itMs)) '_' prms.thrScore{itMs} ]) = thrTmpMap;
    
    elseif ~thrByMapMsInd( itMs )
        % By population measure(s)
        prmsCopy_byPop                         = prms;
        prmsCopy_byPop.thrScore                = prms.thrScore( itMs );
        prmsCopy_byPop.thrPrc                  = prms.thrPrc( itMs );
        [BVCIndByTrial_MsTemp, thrTmpPop]      = classifyByPopSubFunc( Res, wellSampledInd, prmsCopy_byPop );
        rThrAll( itMs, :, : )                  = thrTmpPop;
    end

    % Unify classification for each measure.
    BVCIndByTrial  = BVCIndByTrial & BVCIndByTrial_MsTemp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If requested, remove BVCs better and significantly fit by place cell model. Fields 'isPCByTrial' and 
% 'PCRespMax' need to already be present in table, these are generated by function 'analyseBVCvsPCFits'.
if prms.excPCFits
    bstFitPCInd   = Res.isPCByTrial(:,prms.thrTrials) & Res.PCRespMax(:,prms.thrTrials)>Res.BVCRespMax(:,prms.thrTrials);
    BVCIndByTrial = BVCIndByTrial & ~bstFitPCInd;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert the 'by rate map' classifications to 'by cell'. (And other final clean up).
% Go from two trials to one clasiifcation, using appropriate AND or OR logic %
if strcmp( prms.thrTrialAv, 'min' )
    BVCInd = all( BVCIndByTrial, 2 );
elseif strcmp( prms.thrTrialAv, 'max' )
    BVCInd = any( BVCIndByTrial, 2 );   
end   

% Assign final output %
Res.( prms.outputFieldName )             = BVCInd;
Res.( [prms.outputFieldName 'ByTrial'] ) = BVCIndByTrial;

% Print thresholds
disp(rThrAll(:,:,1)); disp(rThrAll(:,:,2));

% Delete shuf score fields so Res doens't become huge.
for itFd=1:length(ShufScores)
    Res.(ShufScores{itFd}) = [];
end
% Also delete BVCRespMax_VFS field, which was only there to make the naming convention work.
if any(strcmp(ShufScores,'ShufBVCRespMax_VFS'))
    Res.BVCRespMax_VFS = [];
end

%%%%%%%%%%%%%%%%%%%%%%%
end % Main Function END
%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BVCIndByTrial, rThrMed, rThrAll] = classifyByMapSubFunc( Res, wellSampledInd, prms )
%%% Calculate thresholds and classify map-by-map.
tmpClassDiffMs = true(size(Res,1),length(prms.thrTrials));
for itMs=1:length( prms.thrScore ) % itMs= iterator for measure. NOTE: I think this is redundant now, all calls should have only one thrMeasure specified.
    
    % Get shuf score percentiles, real scores, and compare
    shufScores       = Res.( ['Shuf' prms.thrScore{itMs}] )(:,prms.thrTrialsShuf);
    shufScores( cellfun(@isempty, shufScores) ) = {nan};    % Otherwise PRCTILE won't operate via CELLFUN
    rThrAll             = cellfun( @(x,p) prctile(x,prms.thrPrc(itMs)), shufScores);
    rTest            = Res.( prms.thrScore{itMs} )( :, prms.thrTrials );
    tmpClassDiffMs   = tmpClassDiffMs  &  rTest>=rThrAll;
    
    % Keep a record for info of thresholds (median for each age bin group)
    for itAB=1:size(prms.ageBins,1)
        for itBR = 1:2
            rThrForGroup            = rThrAll(  Res.age>=prms.ageBins(itAB,1) & Res.age<=prms.ageBins(itAB,2) & Res.area==itBR,    :    );
            if isempty( rThrForGroup )
                 rThrMed(itMs,itAB,itBR) = nan;
            else
                rThrMed(itMs,itAB,itBR)  = median(rThrForGroup,"all","omitnan");
            end
        end
    end
end
BVCIndByTrial                  = tmpClassDiffMs;
% Insufficently well-sampled cells (on specific trials) cannot be BVCs:
BVCIndByTrial( ~wellSampledInd(:,prms.thrTrials) ) = false;

end % Sub-func END


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BVCIndForAllPopGr, rThrAll, isSpat]  = classifyByPopSubFunc( Res, wellSampledInd, prms )
%%% Classify using a brain area and region-specific set of population thresholds.
BVCIndForAllPopGr = false(size(Res,1),length(prms.thrTrials));
isSpat            = false(size(Res,1),1);
for itAB=1:size(prms.ageBins,1)   
    for itBR = 1:2

        % Pre-allocate temp classification index for this age bin/brain region group (to allow multiple classifcation measures, start with all true, then test with AND in the loop below).
        groupInd         = Res.age>=prms.ageBins(itAB,1) & Res.age<=prms.ageBins(itAB,2) & Res.area==itBR;
        tmpClassDiffMs   = repmat( groupInd, 1, length(prms.thrTrials) );  % Index for combining classification across different measures.
        tmp_isSpat = groupInd;
        
        %%% Loop through classification measures ..
        for itMs=1:length( prms.thrScore ) % itMs= iterator for measure
        
            %%% Get the test scores (i.e. real data) %%%
            rTest                                      = Res.( prms.thrScore{itMs} )( :, prms.thrTrials );
            rTest( ~wellSampledInd(:,prms.thrTrials) ) = nan;   % Cells not sufficently well sampled cannot be BVCs.
        
            %%% Generate thresholds (or if fixed thresholds supplied, get the appropriate value for age and measure) %%%
            if isempty( prms.fixedThrs )

                %%% Normal operation: get thresholds from shuffled data %%%%
                shufDataforThrTrs                  =  Res.( ['Shuf' prms.thrScore{itMs}] )(:,prms.thrTrialsShuf);
                shufDataforThrTrs( ~wellSampledInd(:,prms.thrTrialsShuf) ) = {[]};    % Remove insuffciently sampled cells from random population. 'wellSampledInd' is the same n cols as prms.thrTrials and prms.thrTrialsShuf.
                shufDataforThrTrs                  = shufDataforThrTrs( groupInd, : );
                shufScoresVect                     = cell2mat( reshape( shufDataforThrTrs, [], 1 )  );   % Need to reshape to a column, as sometimes there are empty arrays in the cell array.
                rThr                               = prctile( shufScoresVect(:), prms.thrPrc(itMs) );
                rThrAll(itMs,itAB,itBR)            = rThr;  % Store all thrs for later display.
                
            else
                % Fixed thresholds: just get the relevant value from input param struct %
                rThr = prms.fixedThrs( itMs, itAB, itBR );
            end

            %%% Test scores for measure itMs against the threshold %%%
            tmpClassDiffMs = tmpClassDiffMs  &  rTest>=rThr;
            
            if strcmp(prms.thrScore{itMs},'SI')
                tmp_isSpat = tmp_isSpat & any(rTest>=rThr,2);
            end

            %%% END Loop itMs - different classifier measures %%%
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Combine the classifier index for this age bin/brain region with the full classifier array % 
        BVCIndForAllPopGr = BVCIndForAllPopGr | tmpClassDiffMs;
        
        isSpat = isSpat | tmp_isSpat;

        %%%% END for itBR %%%
    end
    %%%% END for itAB %%%%
end

end % Sub-func end


