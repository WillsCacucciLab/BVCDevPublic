function [ Res ] = classifyCellType( Res, ResShuf, varargin )
% This is a wrapper for running the BVC, BVC (in CCE) and BC classification on whole dataset.
%
% Calls 'shufThrCellClassifer' three types, with parameters set to classify each 
% of the three categories. Also constructs a general filter index, stored as Res.dataInd.

% General params:
prms.classify.BVC.ageBins               = [16 18; 19 21; 22 25; 40 40]; % copy over
prms.classify.BVC.thrTrials             = 1:2;    % Trials whose actual real data scores are used to define BVCs
prms.classify.BVC.thrTrialsShuf         = 1:2;    % Trials from which shuffled data is taken
prms.classify.BVC.thrTrialAv            = 'max'; % 'sortFirst';   'mean'; % 'max';'min;   % How to get an average test value when more than one trial is specified in 'prms.thrTrials'. % default: 'max' (also not important when using 1 trial only for def)
prms.classify.BVC.useOnly2TrBsl         = 0;
prms.classify.BVC.meanRateThr           = 0.2; % default: 0.2 (Bjerknes paper)
prms.classify.BVC.nBinVisThr            = 0.8 .* [625 784 804];  % 0.8 [627 specific for Moser data: 0.8 * (28^2)]
prms.classify.BVC.filterByWaveformProps = 0; % default: no filter
prms.classify.BVC.minSpkWidth           = 0.3; % in ms - waveforms have to be wider than this (peak-to-trough)
% BVC specific:
prms.classify.BVC.thrByMap             = [0 0];    % If 1, then thresholds are calculated for individual maps, and BVCs defined on if their scores crosses their self-generated threshold 
prms.classify.BVC.thrScore             = {'BVCRespMax','SI'}; % 'BVCRespMax'; %   % default: {}
prms.classify.BVC.thrPrc               = [99 75]; % 99; %  99; %  default: 99.75 (in combination with next 2 param defaults)
prms.classify.BVC.excPCFits            = 0; %% Leaving this here as we might need for upcoming analysis.
prms.classify.BVC.fixedThrs            = [];
% Border cell specific:
prms.classify.BC                   = prms.classify.BVC;   % Copy over to get copy of general params
prms.classify.BVC.thrByMap         = [0 0];
prms.classify.BC.thrScore          = {'borderScore','SI'}; % default: {'borderScoreAdSm_New','SI'} (Bjerknes paper)
prms.classify.BC.thrPrc            = [95 95]; % default: 95 (Bjerknes paper)
prms.classify.BC.loadBCsFromXCL    = 0;
prms.classify.BC.path2TalesList    = 'D:\Dropbox\2017_02 Sub Dev project\';

% -----------------------------------------------------------------------------------%
% - Because of the complex nature of the params struct, here you either supply the   %
%   struct, or you call this function and use the parameters written above: you      %
%   can't give single options as input arg pairs.                                    %
if ~isempty(varargin)                                                                %
    prms = varargin{1};                                                               %
end                                                                                  %
% ---------------------------------------------------------------------------------- %

% If requested, filter by spikewidth %
if prms.BVC.filterByWaveformProps
    % Filter out narrow waveforms (put. interneurons)
    spkWdthInd = Res.spkWidthMean >= prms.minSpkWidth;
    Res = Res(spkWdthInd, : );
    ResShuf = ResShuf(ismember(ResShuf.cellID,Res.cellID), : ); % different oder in table compared to 'Res'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classify cells as BVC, BC, VRC.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) BVCs (i.e. using classic brute force search of best fit model).
fprintf('\n====== Thrs for BVCs\n'); %: %s prc=%4.2f, SUB & mEC ==========\n', prms.thrScoreBVC, prms.thrPrcBVC);
prms_cp                 = prms.BVC;
prms_cp.outputFieldName = 'isBVC';
[Res, rThr]             = shufThrCellClassifier( Res, ResShuf, prms_cp );
% Save classification params.
Res.Properties.UserData.currClass.BVC          = prms.BVC;
Res.Properties.UserData.currClass.BVC.rThr     = rThr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) BVCs, defined in the CCE trial.
fprintf('\n====== Thrs for BVCs in CCE:\n'); % %s prc=%4.2f, SUB & mEC ==========\n', prms.thrScoreBVC, prms.thrPrcBVC);
prms_cp.outputFieldName = 'isBVC_CCE';
prms_cp.thrTrials       = 5;
prms_cp.thrTrialsShuf   = 5;
[Res, rThr]             = shufThrCellClassifier( Res, ResShuf, prms_cp );
% Save classification params.
Res.Properties.UserData.currClass.BVC_CCE          = prms.BVC;
Res.Properties.UserData.currClass.BVC_CCE.rThr     = rThr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Border cells. i.e. normally SI and border score. For Moser data can load definition from list.
if prms.BC.loadBCsFromXCL 
    Res = loadFromXCL(Res,0,prms.BC );
    fprintf('%s\n','================ BC data loaded from Tale''s list ================');
else
    fprintf('%s\n','=================== Thresholds for BCs in SUB & mEC ===================\n');
    prms_cp                 = prms.BC;
    prms_cp.outputFieldName = 'isBC';
    [Res, rThr]             = shufThrCellClassifier( Res, ResShuf, prms_cp );
end
Res.Properties.UserData.currClass.BC          = prms.BC;
Res.Properties.UserData.currClass.BC.rThr     = rThr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Spatial cells (i.e. just those crossing a SI threshold.
fprintf('\n====== Thrs for Spat\n'); %: %s prc=%4.2f, SUB & mEC ==========\n', prms.thrScoreBVC, prms.thrPrcBVC);
prms_cp                 = prms.spat;
prms_cp.outputFieldName = 'isSpat';
[Res, rThr]             = shufThrCellClassifier( Res, ResShuf, prms_cp );


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (4) HD cells.
% fprintf('\n====== Thrs for HDCs\n'); %: %s prc=%4.2f, SUB & mEC ==========\n', prms.thrScoreBVC, prms.thrPrcBVC);
% prms_cp                 = prms.HD;
% prms_cp.outputFieldName = 'isHD';
% [Res, rThr]             = shufThrCellClassifier( Res, ResShuf, prms_cp );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a general filter index, that gets stored in Res.dataInd.
% NB TW - I'm not sure what this index gets used for, but I am preserving it in case important.
% Make the a general index for all the data
nBinVisThrAbs                = nan(height(Res),1);
nBinVisThrAbs(Res.area == 1) = prms.BVC.nBinVisThr(1);
nBinVisThrAbs(Res.area == 2) = prms.BVC.nBinVisThr(2);
sampleFilterForAll(:,1:4)    = Res.meanRate(:,1:4) >= prms.BVC.meanRateThr & Res.nBinVis(:,1:4) >= nBinVisThrAbs;
sampleFilterForAll(:,5)      = Res.meanRate(:,5) >= prms.BVC.meanRateThr & Res.nBinVis(:,5) >= prms.BVC.nBinVisThr(3);
% TW: I added in the filter for cross-tet recorded cells here, it seems like a sensible place. 
% The index 'isXRecd' is generated by function 'crossRecCellsChecker'.  
if ~ismember('isXRecd', Res.Properties.VariableNames); Res = crossRecCellsChecker(Res); end
sampleFilterForAll(Res.isXRecd, :) = false;
% If requested, remove data with only one BSL trial (only relevant to EC data) %
if prms.BVC.useOnly2TrBsl
   addInd = ~isnan( Res.trialInd(:,2));
else
   addInd = true(height(Res),1);
end
% Combine filters and save.
Res.dataInd                  = sampleFilterForAll & addInd;


