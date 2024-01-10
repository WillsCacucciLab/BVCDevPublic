function [FigPl, FigHst, statsOut] = Figure3HistogramsAndSummaryPlots(Res,varargin)
% New function for MS, just produces plots and stats for figure 3, nothing extra whatsoever.


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

FigPl   = gra_multiplot( length(prms.measuresToPlot), 1,                    'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25, 'structOut', 1 );   axArrPl=FigPl.axArr;
FigHst  = gra_multiplot( length(prms.measuresToPlot), size(prms.ageBins,1), 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25, 'structOut', 1 );   axArrHst=FigHst.axArr;
% hFigScat = gra_multiplot( 1,                           size(prms.ageBins,1), 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25 );   axArrScat=getappdata(hFigScat,'axesHandles');

dAndPhiForPy = cell( 2, size(prms.ageBins,1) );

for itMs=1:length(prms.measuresToPlot)

    [M,E]         = deal(nan(size(prms.ageBins,1),1));
    E             = repmat(E, [1 1 2]);  % E needs 3rd dim to store asymmetric errors for barwitherr 
    dataForStats  = [];
    grVarForStats = ones(0,2);
    
    for itAB = 1:size(prms.ageBins,1)
        for itSp=1:size(prms.dSplit,1)

            % Get data for age and cell type group. Selection of data trial-by-trial can depend on measure.
            indForGr  = Res.age >= prms.ageBins(itAB,1) & Res.age <= prms.ageBins(itAB,2) & Res.area==prms.area & Res.(prms.cellType) & any( Res.dataInd(:,2), 2 ); 
            % Need to select single 'best' BVC baseline trial.
            dataForGr     = Res.( prms.measuresToPlot{itMs} )( indForGr, prms.bslTrToPlot );
            if strcmp(prms.cellType,'isBVC')
                [~,bstBslCol] = max( Res.BVCRespMax( indForGr, prms.bslTrToPlot ), [], 2, "omitnan" );
            elseif strcmp(prms.cellType,'isBC')
                [~,bstBslCol] = max( Res.borderScore( indForGr, prms.bslTrToPlot ), [], 2, "omitnan" );
            end
            bstBslLin     = sub2ind(size(dataForGr), 1:size(dataForGr,1), bstBslCol');
            dataForGr     = dataForGr(bstBslLin');

            % Always need to get dMax, so as to split other measures by dMax.
            dMaxForGr     = Res.dMax( indForGr, prms.bslTrToPlot );
            dMaxForGr     = dMaxForGr(bstBslLin');
            dMaxSplitInd  = dMaxForGr >= prms.dSplit(itSp,1) & dMaxForGr <= prms.dSplit(itSp,2);
            dataForGr     = dataForGr( dMaxSplitInd );

            % Plot histogram for this age bin.
            if contains(prms.measuresToPlot{itMs},'dMax')
                dataForGr = ( dataForGr-1 ) .* 2.5;
                histBins  = 0:2.5:32.5;
            elseif contains(prms.measuresToPlot{itMs},'phiMax') && ~prms.phiHistQuadAvs
                dataForGr = mod(dataForGr,pi*2);
                histBins  = linspace(0,2*pi,61);
                crnrLines = (0:(pi/2):(3*pi/2)) + pi/4;  % For later plotting
            elseif contains(prms.measuresToPlot{itMs},'phiMax') && prms.phiHistQuadAvs
                dataForGr = mod(dataForGr+pi/4,pi/2);
                histBins  = linspace(0, pi/2, 16);
                crnrLines = pi/4;
            elseif contains(prms.measuresToPlot{itMs},'SZMaxInd')
                if prms.plotSZ_HHW
                    dataForGr = transpose( Res.Properties.UserData.sz_hhw_list( dataForGr' ) );
                    histBins  = 0:0.5:15;
                else
                    histBins  = 0.5:1:4.5;
                end
            end
            ax = axArrHst( itMs, itAB );
            histogram( ax, dataForGr, histBins, 'DisplayStyle','bar','edgecolor','none' );
            % Histogram formatting
            if contains(prms.measuresToPlot{itMs},'phiMax')
                hold(ax,'on');
                for ii=crnrLines;  plot(ax, [1 1].*ii, ax.YLim, 'k:');   end
                if ~prms.phiHistQuadAvs
                    tickPos = (0:(pi/2):(3*pi/2))+pi/4;
                    tickLbl = {'0','90','180','270'};
                else
                    tickPos = [0 pi/4 pi/2];
                    tickLbl = {'-45','0','45'};
                end
                set(ax,'xTick', tickPos, 'xticklabel', tickLbl, 'xticklabelrotation',45);
            end
            ylabel(ax,'Cell count');
            xlabel(ax, prms.xLabels{itMs}, 'Interpreter','tex');

            % Scatter plot of phi vs d.
            % - In the end I needed to do this plot in python, so all we do here is package data for export.
            if contains(prms.measuresToPlot{itMs},'phiMax')
                dForPlot = 14-dMaxForGr;
                dAndPhiForPy{1,itAB} = dForPlot;
                dAndPhiForPy{2,itAB} = dataForGr;
            end

            % Get group data average and spread measures
            if any( contains({'dMax','SZMaxInd'}, prms.measuresToPlot{itMs}) )
                % Max D of BVC tuning obviously non-normal, also SZ just odd, so use median +-quartiles.
                % BARWITHERR will do asymmetric error bars, stack in 3rd dim.
                M(itAB,itSp)    = median(dataForGr,1,'omitnan');
                E(itAB,itSp,1)  = prctile(dataForGr,25);
                E(itAB,itSp,2)  = prctile(dataForGr,75);
                E(itAB,itSp,:)  = E(itAB,itSp,:) - M(itAB,itSp);  % Convert absolute error limits to err bar lengths.
            elseif contains(prms.measuresToPlot{itMs},'phiMax')
                % Phi Max (Part 1): test RV significance after angle quadrupling.
                phiQuad         = mod( dataForGr.*4, pi*2 );     % Data has been shifted by 45deg, but not relevant to this test.
                [p,z]           = circ_rtest(phiQuad);
                fprintf(1,'\nAge bin %d R-test (quad), z=%3.1f, p=%4.3f', itAB, z, p);
                statsOut.rQuad_z(itAB) = z;
                statsOut.rQuad_p(itAB) = p;
                % Phi Max (Part 2): measure for development summary plot is proportion of phi's within +-12deg range of wall ori.
                dataForGr       = mod( dataForGr, pi/2);                      % Already shifted 0deg -> 45deg, for histogram above.
                wallOriInd      = dataForGr>32/180*pi & dataForGr<58/180*pi;  % This is the 45deg bin and the two flanking on each side (+-12deg). Go a degree further than this, to avoid missing data through rounding issues. 
                propWall        = sum(wallOriInd) / length(dataForGr);
                M(itAB,itSp)    = propWall.*100;   
                E(itAB,itSp,:)  = 1.96  .*   sqrt(   ( propWall.*(1-propWall) ) ./ length(dataForGr)  )   ./   2   .*100;  % Zar Eq 24.63
                dataForGr       = double(wallOriInd); % Turn data into binomial wall/no-wall counts: stats test for this is chi-sq (below)
            end
            dataForStats  = [dataForStats; dataForGr];
            grVarForStats = [grVarForStats;  [ones(length(dataForGr),1).*itAB, ones(length(dataForGr),1).*itSp] ];
        
        end
    end
    
    % Summary plots, i.e. average and error for each age bin.
    ax            = axArrPl(itMs);  axes(ax);
    [hBar, ~]     = barwitherr(E,M);
%     hBar.BarWidth = deal( 0.7 );
    set(ax,'xticklabel',prms.ageLabels,'xticklabelrotation',45,'box','off');
    ylabel(ax,prms.yLabels{itMs}, 'Interpreter','tex');
    if contains(prms.measuresToPlot{itMs},'phiMax')
        hold(ax,'on');
        plot(ax.XLim, [33.33 33.33],'k:');
    end

    % Stats. Run one 1-way ANOVA (or equivalent) per measure. 
    fprintf(1,'\n%s\n',prms.measuresToPlot{itMs});
    if contains({'dMax'}, prms.measuresToPlot{itMs})
        % dMax needs KW test as non-normal. Output 'T' same format as ANOVA1, except chi-sq takes place of F.
        [~, T]        = kruskalwallis(dataForStats, grVarForStats(:,1), 'off');
    elseif contains({'SZMaxInd'}, prms.measuresToPlot{itMs})
        % SZMaxInd weird - currently using KW when just age, but chi-sq when age*d - just because I don't know how to do the latter better.
        if size(prms.dSplit,1)==1
            % Regular 'all data' analysis (no split by d tuning)
            [~, T]        = kruskalwallis(dataForStats, grVarForStats(:,1), 'off');
        elseif size(prms.dSplit,1)==2
            % This is when data is split into long and short d tunings - add d tuning category to crosstab analysis.           
            for itAB=1:4
                ABInd = grVarForStats(:,1)==itAB;
                P     = ranksum( dataForStats(ABInd), grVarForStats(ABInd,2) );
                disp(P);
            end 
            % WHEN RUNNING THIS, DON'T CONVERT SZ INDICES TO HHWs!
            [~,chiSq,p] = crosstab(dataForStats, grVarForStats(:,1), grVarForStats(:,2));
            % Fake an ANOVA output for printing compatibility
            T = cell(3,6);  T{2,6}=p;  T{2,5}=chiSq;  T{2,3}=0; T{3,3}=0;
        end
    elseif contains(prms.measuresToPlot{itMs},'phiMax')
        % phiMax is actually, by this point, binomial counts of wall vs non-wall oriented. Use Age*(Wall/non-wall) Chi-sq.        
        if size(prms.dSplit,1)==1
            % Regular 'all data' analysis (no split by d tuning)
            [~,chiSq,p] = crosstab(dataForStats, grVarForStats(:,1));
        elseif size(prms.dSplit,1)==2
            % This is when data is split into long and short d tunings - add d tuning category to crosstab analysis.
            [tbl,chiSq,p] = crosstab(dataForStats, grVarForStats(:,1), grVarForStats(:,2));
            statsOut.phiSplitChi = chiSq;
            statsOut.phiSplit_p  = p;
            % 'Post-hoc' z-test comparison of proportions, with AB
            P = transpose( tbl(:,:,1) );
            N = transpose( sum(tbl,3) );
            fprintf( 1, '\nZ-test proportions, wall vs non-wall within Age bin:\n' );
            for itZT=1:size(P,1)
                [~,p_zt, chi2stat_zt,df_zt] = prop_test(  P(itZT,:),  N(itZT,:), 0  );
                fprintf(1, 'Age Bin=%d, Chi2(%d)=%3.2f, p=%4.3f\n', itZT, df_zt, chi2stat_zt, p_zt );
                if itZT==1
                    statsOut.phiSplitP16_z  = chi2stat_zt;
                    statsOut.phiSplitP16_df = df_zt;
                    statsOut.phiSplitP16_p  = p_zt;
                end
            end

        end
        % Fake an ANOVA output for printing compatibility
        T = cell(3,6);  T{2,6}=p;  T{2,5}=chiSq;  T{2,3}=0; T{3,3}=0;
    end
    fprintf(1,'D-range %d-%d: F(%d,%d)=%3.2f, p=%4.3f; \t', prms.dSplit(itSp,1), prms.dSplit(itSp,1), T{2,3}, T{3,3}, T{2,5}, T{2,6});

end

save( 'dAndPhiForPy.mat', 'dAndPhiForPy', '-v7');  % Save d and phi data (for 2D rad hists). Make plots in python.
gra_multilabel(FigHst.hFig,'col',prms.ageLabels);






