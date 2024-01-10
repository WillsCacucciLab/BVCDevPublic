function Res = finalFigurePlots_rev1( Res, ResShuf, modelMaps, varargin )
% Classify BVCs and plot basic spatial props %

% Load global parameter file
prms = globalParams_final;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - which figures to plot  %
prms.plotClassFig                  = 0;  % Figure 1 Panel with histograms of BVC response
prms.calcProportions               = 0;  % Figure 2 - proportions of BVCs
prms.plotCellProps                 = 0;  % Figure 2 - BVC properties by age.
prms.figure3Hists                  = 0;
prms.figure4Plots                  = 1;
prms.doBarrierAnalysis             = 0;  % Final MS figure 5 plots are within this figure.
prms.BVCvsBCAnalysis               = 0;  % Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Classify the whole dataset
[ Res ] = classifyCellType( Res, ResShuf, prms.classify );    


%% Plot Classification scores figure
% Plots of BVC and vector map scores, BVCs in blue non-BVC in black
if prms.plotClassFig
    plotClassificationFigure( Res, ResShuf, modelMaps );
end

%% Plot Cell Type Proportions
% plot proportions for cell type in selected area (and plot)
if prms.calcProportions
    getCellTypeProportions( Res, prms.area,'ageBins', prms.ageBins, 'ageLabels', prms.ageLabels );
end

%% Plot Some Basic Stats across age
if prms.plotCellProps
    if ~isempty(prms.fig2props.subSampMode)
        Res = subSampleRes( Res, prms.fig2props );  % Only used ATM for excluding drifty clusters ('wfPkDiff')
    end
    plotCellPropsByAge( Res, prms.area, prms.fig2props );
end

%% Figure 3 - histograms and analysis of BVCs in square
if prms.figure3Hists
    Figure3HistogramsAndSummaryPlots( Res, prms.fig3 );
end

%% Figure 4 - distributions of phi in CCE
if prms.figure4Plots
    figure4Plots(Res,prms.fig4);
end


%% Figure 5 - barrier trials
if prms.doBarrierAnalysis
    prms.barrAnalysis.modMaps = modelMaps;
    prms.barrAnalysis.area    = prms.area;   
    barrierAnalysis_v2( Res, prms.area, prms.barrAnalysis );
end


%% Comparison of BVCs vs BCs
if prms.BVCvsBCAnalysis
    compareBVCvsBC(Res,prms.BVCvsBC);
end






