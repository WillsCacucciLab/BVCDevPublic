function [] = compareBVCvsBC(Res,varargin)
% Plot and stats to compare development of BVCs and BCs in a given brain region.
% Wrote this as 'plotCellPropsByAge' was getting to long and unweildy to edit easily,
% this function should be shorter and more focused.

% TODO, 01/08/22
% Change AR to be longer, for more age bins
% Line plots not bar plots for first three plots
% Assymetric error bars for d

% Params
prms.ageBins   = [16 18; 19 21; 22 24; 25 27; 28 30; 31 33; 34 36; 40 40];  % [16 18; 19 21; 22 25; 26 36; 40 40];  %  default: [16 18; 19 21; 22 25; 26 36; 40 40];     [15 18; 19 21; 22 24; 25 27; 28 30; 31 33; 34 36; 40 40];  %   [15 17; 18 19; 20 21; 22 25; 26 36; 40 40];   %    
prms.ageLabels = {'P16-18','P19-21','P22-24','P25-27','P28-30','P31-33','P34-36','Adult'};
prms.area                 = 2;
prms.cellTypeList         = {'isBVC','isBC'};
prms.bslTrToPlot          = 1:2;
prms.plotOnlyIsBVCTrials  = 0;

prms.measuresToPlot = {'SI',   'BVCRespMax', 'interStab','intraStab',  'dMax', 'phiMax'}; %
prms.yLimMax        = {[0 1.5], [0.4 1],      [-0.2 1] ,  [-0.2 1],    [0 20], [0 100]}; %,   
prms.yLabels        = { 'Spatial Information (bits/spk)','BVC r_m_a_x', 'Inter-trial Stability (r)','Intra-trial stability (r)','{\itd} (cm)',   ['% Wall oriented ' char(981)]  };
lSpecForCT          = {'b-','r-'};
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

hFigPl  = gra_multiplot( length(prms.measuresToPlot), 1, 'plotsize', [4 2], 'axesborder', [1 1 1 1].*1.25 );   
axArr  = getappdata(hFigPl,'axesHandles');
hFigHst = gra_multiplot( 1, size(prms.ageBins,1), 'plotsize', [4 2], 'axesborder', [1 1 1 1].*1.25 );   
axArrHst  = getappdata(hFigHst,'axesHandles');


for itMs=1:length(prms.measuresToPlot)

    [M,E]         = deal(nan(size(prms.ageBins,1),2));
    dataForStats  = [];
    grVarForStats = ones(0,2);
    
    for itAB = 1:size(prms.ageBins,1)
        for itCT = 1:2
            
            % Get data for age and cell type group. Selection of data trial-by-trial can depend on measure.
            indForGr  = Res.age >= prms.ageBins(itAB,1) & Res.age <= prms.ageBins(itAB,2) & Res.area==prms.area & Res.(prms.cellTypeList{itCT}); 
            cellIDForGr = Res.cellID( indForGr );
            SIForGr     = Res.SI( indForGr );
            if contains(prms.measuresToPlot{itMs},'inter')
                % Inter-trial stabilities - only one inter-bsl comparison possible, so just one column.
                dataForGr = Res.( prms.measuresToPlot{itMs} )( indForGr, 1 );
            elseif contains(prms.measuresToPlot{itMs},'phiMax')
                % Phi max - need to select single 'best' BVC trial 
                dataForGr     = Res.( prms.measuresToPlot{itMs} )( indForGr, prms.bslTrToPlot );
                [~,bstBslCol] = max( Res.BVCRespMax( indForGr, prms.bslTrToPlot ), [], 2, "omitnan" );
                bstBslLin     = sub2ind(size(dataForGr), 1:size(dataForGr,1), bstBslCol');
                dataForGr     = dataForGr(bstBslLin');
            else
                % This is the normal case for most measures, with the exception that ..
                % .. for border cells, BVCRespMax should be replaced with border score, so the two
                % different measures can be on one plot.
                if contains(prms.measuresToPlot{itMs},'BVCResp') && strcmp( prms.cellTypeList{itCT}, 'isBC' )
                    dataForGr = Res.borderScore( indForGr, prms.bslTrToPlot );
                else
                    % This is actually the real default, normal case.
                    dataForGr = Res.( prms.measuresToPlot{itMs} )( indForGr, prms.bslTrToPlot );
                end
                if prms.plotOnlyIsBVCTrials
                    trialSelInd             = Res.([prms.cellTypeList{itCT} 'ByTrial'])( indForGr, prms.bslTrToPlot );
                    dataForGr(~trialSelInd) = nan; 
                end
                % Select 'best' baseline trial, for all measures.
                if strcmp(prms.cellTypeList{itCT},'isBVC')
                    [~,bstBslCol] = max( Res.BVCRespMax( indForGr, prms.bslTrToPlot ), [], 2, "omitnan" );
                elseif strcmp(prms.cellTypeList{itCT},'isBC')
                    [~,bstBslCol] = max( Res.borderScore( indForGr, prms.bslTrToPlot ), [], 2, "omitnan" );
                end
                bstBslLin     = sub2ind(size(dataForGr), 1:size(dataForGr,1), bstBslCol');
                dataForGr     = dataForGr(bstBslLin');
            end
            % Make sure dMax index is converted to actual cm!
            if strcmp(prms.measuresToPlot{itMs},'dMax')
                dataForGr = (dataForGr-1) .* 2.5;
            end
            
            if itAB==8 && strcmp(prms.measuresToPlot{itMs},'dMax') % 
                ind = find(dataForGr>=10);
                for ii=1:length( ind )
                    fprintf(1, '\n %s \t %4.3f \t %2.1f', cellIDForGr{ind(ii)}, SIForGr(ind(ii)), dataForGr(ind(ii)));
                end
            end

            % Histograms for BVC d for Edvard
            if strcmp(prms.cellTypeList{itCT},'isBVC') && strcmp(prms.measuresToPlot{itMs},'dMax')
                histBins  = 0:2.5:32.5;
                histogram( axArrHst(itAB), dataForGr, histBins,  'DisplayStyle','bar','edgecolor','none' );
                ylabel(axArrHst(itAB),'Cell count');
                xlabel(axArrHst(itAB), '{\itd} (cm)', 'Interpreter','tex');
            end
            
            % Get group data average and spread measures
            if contains(prms.measuresToPlot{itMs},'dMax')
                % Max D of BVC tuning obviously non-normal needs use median +-IQR
                % Note on error bars - really need a box-plot here, but for convenience calculating
                % the upper quartile, which is the useful measure of spread. In all cases except
                % one BVC age bin, 25th %ile = 1. Need to change function to get final plot.
                M(itAB,itCT)   = median(dataForGr,1,'omitnan');
                E(itAB,itCT,1) = prctile(dataForGr,50) - prctile(dataForGr,25);
                E(itAB,itCT,2) = prctile(dataForGr,75) - prctile(dataForGr,50);
            elseif contains(prms.measuresToPlot{itMs},'phiMax')
                % Phi Max - measure calculated is proportion of phi's within +-9deg range of wall ori.
                dataForGr      = mod( dataForGr+pi/4, pi/2);
                wallOriInd     = sum(dataForGr>=32/180*pi & dataForGr<58/180*pi);
                propWall       = wallOriInd / length(dataForGr);
                M(itAB,itCT)   = propWall*100;   
                E(itAB,itCT)   = 1.96  .*   sqrt(   ( propWall.*(1-propWall) ) ./ length(dataForGr)  )   ./   2  .*  100;  % Zar Eq 24.63
                if E(itAB,itCT)<=100;   E(itAB,itCT,2) = E(itAB,itCT,1);   else;    E(itAB,itCT,2) = 100;   end
                dataForGr      = double(wallOriInd); % Turn data into binomial wall/no-wall counts: stats test for this is chi-sq (below)
            else
                % This is the default case:
                M(itAB,itCT)  = mean(dataForGr,1,'omitnan');
                E(itAB,itCT)  = std(dataForGr,'omitnan') / sqrt( sum(~isnan(dataForGr)));
            end
            dataForStats  = [dataForStats; dataForGr];
            grVarForStats = [grVarForStats;  [ones(length(dataForGr),1).*itAB, ones(length(dataForGr),1).*itCT] ];

        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final Manuscript Plots.
    ax = axArr(itMs);
    if contains(prms.measuresToPlot{itMs},'Phi') 
        axes(ax);
        [hBar, ~]     = barwitherr(E,M);
        set(ax,'xticklabel',prms.ageLabels);
    else
        if 0
            xOffset = [0 0];
            for itCT=1:2
                errorbar(ax, (1:size(M,1))+xOffset(itCT), M(:,itCT), E(:,itCT), lSpecForCT{itCT});
                hold(ax,'on');
            end
            set(ax, 'xlim', [0.5 size(M,1)+0.5], 'xtick', 1:size(M,1));
        else
            d4Box            = dataForStats;
            gr4Box           = grVarForStats;
            labels           = [prms.ageLabels; repmat({''},1,length(prms.ageLabels))];
            labels           = transpose( labels(:) );
            hBP              = boxplot(ax,  d4Box, gr4Box, 'Notch','off','Colors', [0 0.1 1; 1 0 0], 'Symbol', '',...
                               'FactorGap', [3 0], 'boxstyle', 'outline', 'width',0.7,'LabelVerbosity','majorminor'...
                                   , 'plotstyle', 'compact', 'labels', labels);
            set(hBP, 'LineWidth', 0.8);
            hold( ax, 'on' );
            if prms.yLimMax{itMs}(1) < 0
                plot(ax, ax.XLim, [ 0 0 ], 'k:' );
            end
        end
    end
    set(ax,'xticklabelrotation',45,'ylim',prms.yLimMax{itMs},'box','off');
    ylabel(ax,prms.yLabels{itMs}, 'Interpreter','tex');
    if contains(prms.measuresToPlot{itMs},'phiMax')
        hold(ax,'on');
        plot(ax.XLim, [33.33 33.33],'k:');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Stats: test for age effects.
    % Originally wrote this do put all data into a 2-way ANOVA, but can't actually do that as cells that are both 
    % BVCs and BCs are resampled. Run one 1-way ANOVA per cell type, therefore. 
    fprintf(1,'\n%s\n',prms.measuresToPlot{itMs});
    for itCT=1:2
        indForCT = grVarForStats(:,2)==itCT;
        if contains(prms.measuresToPlot{itMs},'dMax')
            % dMax needs KW test as utterly non-normal. Output 'T' same format as ANOVA1, except chi-sq takes place of F.
            [~, T]      = kruskalwallis(dataForStats(indForCT), grVarForStats(indForCT,1), 'off');
        elseif contains(prms.measuresToPlot{itMs},'phiMax')
            % phiMax is actually, by this point, binomial counts of wall vs non-wall oriented. Use Age*(Wall/non-wall) Chi-sq.
            [~,chiSq,p] = crosstab(dataForStats(indForCT), grVarForStats(indForCT,1));
            % Fake an ANOVA output for printing compatibility
            T = cell(3,6);  T{2,6}=p;  T{2,5}=chiSq;  T{2,3}=0; T{3,3}=0;
        else
            % Otherwise, default test is 1-way ANOVA
            [~, T]      = anova1(dataForStats(indForCT), grVarForStats(indForCT,1), 'off');
        end
        fprintf(1,'%s: F(%d,%d)=%3.2f, p=%4.3f; \t', prms.cellTypeList{itCT}, T{2,3}, T{3,3}, T{2,5}, T{2,6});
    end
    % Stats: test overall mean d for BVCs vs BCs
    if contains(prms.measuresToPlot{itMs},'dMax')
        p = ranksum( dataForStats(grVarForStats(:,2)==1), dataForStats(grVarForStats(:,2)==2) );
        fprintf( 1, '\n Test d BVC vs BC All Ages: p=%4.3f', p);
    end
end

% gra_multilabel(hFigHst,'col',prms.ageLabels);


%%% Venn diagrams of prevalence, overlap %%%
hFig         = gra_multiplot( 1, size(prms.ageBins,1), 'axesborder', [1 1 1 1].*0.25 );   axArray=getappdata(hFig,'axesHandles');
cntsAllAges = nan( 4, size(prms.ageBins,1) );
for itAB = 1:size(prms.ageBins,1)
    % Get proportions
    indForGr  = Res.age >= prms.ageBins(itAB,1) & Res.age <= prms.ageBins(itAB,2) & Res.area==prms.area; 
    A(1)      = sum(Res.isBVC(indForGr));
    A(2)      = sum(Res.isBC(indForGr));
    I         = sum(Res.isBVC(indForGr) & Res.isBC(indForGr));
    cntsAllAges( :, itAB ) = [A'; I; sum(indForGr)];  % Save all props together for stats (as raw counts).
    A         = A./sum(indForGr);
    I         = I./sum(indForGr);
    % Plot venn diagram.
    ax        = axArray(itAB);  axes(ax);
    [H,S]     = venn(A,I);
    set( H(1), 'FaceColor',[0 0 1] );
    set( H(2), 'FaceColor',[1 0 0] );
    axis("equal");
    % Scale these to match absolute percentages to circle sizes across age bins (area inside axis box =100% of cells)
    circ2BoxR = S.CircleArea(1) ./ (diff(ax.YLim) * diff(ax.XLim));
    sf        = sqrt( circ2BoxR / A(1)  );
    ax.XLim   = ax.XLim.*sf;  % NOTE - you only need to change one axis - as axis(equal) has been applied, other axis will automatically scale with it.
    set(ax,'xtick',[],'ytick',[],'box','on');
    % Text labels
    yPos      = ax.YLim(2) - diff(ax.YLim)/10;
    xPos      = ax.XLim(2) - [diff(ax.YLim)*0.2 diff(ax.YLim)*0.5 diff(ax.YLim)*0.8];
    txt       = [A(2) I A(1)] .*100;
    clr       = {[1 0 0],[1 0 1],[0 0 1]};
    for itLb=1:3
        text(xPos(itLb),yPos,[num2str(txt(itLb),'%3.0f') '%'],'parent',ax,'HorizontalAlignment','center','Color',clr{itLb});
    end
end
gra_multilabel( hFig, 'col', prms.ageLabels);
% Stats - z-test of mean pup proportion vs adult proportion.
fprintf(1,'\n');
for itCT=1:3
    nBVC  = [ sum( cntsAllAges(itCT,1:end-1) ),  cntsAllAges(itCT,end) ];
    nCell = [ sum( cntsAllAges(4,1:end-1) ),     cntsAllAges(4,end) ];
    pHat  = nBVC ./ nCell;
    pBar  = sum( nBVC, 2 )  ./ sum( nCell, 2 );   % Eq 24.50
    qBar  = 1 - pBar;
    Z     = (pHat( 1 ) - pHat( 2 ))    /     sqrt(    ((pBar*qBar)/nCell( 1 ))     +     ((pBar*qBar)/nCell( 2 ))    );    % Eq 24.49
    pVal  = (1 - normcdf(abs(Z),0,1)).*2; % times two of the raw cdf is two-tailed probability
    fprintf(1,'Z-test mean pup vs Ad, CT=%d: Z=%3.2f, p=%4.3f \n', itCT, Z, pVal);
end






