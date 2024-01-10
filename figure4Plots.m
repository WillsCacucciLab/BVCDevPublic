function [FigHst, FigRQ, statsOut] = figure4Plots( Res, varargin )
% Analyse distribution of Phi in CCE.


prms.egEnsembles       = {}; % {'r345_P20','r315_P21'}; % {'r372_1406018', 'r372_140609'}; %  {'r345_P16','r315_P18','r345_P20','r156_P23'};  
                       % Ensemble notes: Adult, r372_140609, t7c7, AC 30, tr 2
                       %                                     t1c5, AC 72, tr 1
                       %                        r372_1406018, t3c2, AC 18, tr2
                       %                                      t6c4, AC 42, tr1


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

% Preprare Res table - use only Sub data, and sort by age for ease of plotting later.
Res = Res(Res.area == 1, : );
Res = sortrows(Res, {'age','dataset'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1 - Just plot tuning orientations histograms of BVCs in CCE.
FigHst = gra_multiplot( length(prms.measuresToPlot), 4, 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25, 'structOut', 1 );   axArrHst=FigHst.axArr;
% if prms.phiHistQuadAvs;    hFigQdShHist = gra_multiplot( 1, 4, 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25 );   axArrQdShHist=getappdata(hFigQdShHist,'axesHandles');   end

dAndPhiForPy = cell( 2, size(prms.ageBins,1) );

for itMs = 1:length(prms.measuresToPlot)
    for itAB = 1:size(prms.ageBins,1)
        
        % Get data for age and cell type group. Selection of data trial-by-trial can depend on measure.
        indForGr  = Res.age >= prms.ageBins(itAB,1) & Res.age <= prms.ageBins(itAB,2) & Res.isBVC_CCE  & Res.dataInd(:,5);
        if prms.exclBadCCEBVCs
            indForGr = indForGr   &   Res.BVCRespMax(:,5) > prms.exclBadCCEBVCThrs(itAB);
        end
        dataForGr = Res.( prms.measuresToPlot{itMs} )( indForGr, 5 );
        dMaxForGr = Res.dMax( indForGr, 5 ) - 1;

        % Stats testing:
        % Standard (for this study) quadrupled angle Rayleigh Test
        phiQuad          = mod( dataForGr.*4, pi*2 ); 
        [p,z]            = circ_rtest(phiQuad);
        statsOut.rQuad_z(itAB) = z;
        statsOut.rQuad_p(itAB) = p;
        fprintf(1,'\nAge bin %d R-test (quad), z=%3.1f, p=%4.3f', itAB, z, p);
        % 'Normal' Rayleigh test for unimodal departure from uniformity. 
        [p,z]           = circ_rtest(dataForGr);
        fprintf(1,'\nAge bin %d R-test (unimodal), z=%3.1f, p=%4.3f', itAB, z, p);
%         % Test corner/non-corner ratio at all possible 'corner' angles (anoher test for quad symmetry shifted).
%         if prms.phiHistQuadAvs
%             shiftInd          = rem( (0:14) + (0:14)', 15 ) + 1;
%             phiCntsSh         = hHist.BinCounts(shiftInd);
%             crnRatSh          = sum(phiCntsSh(:,6:10),2)  ./  sum(phiCntsSh,2);
%             [maxCrnRat,bstSh] = max(crnRatSh);
%             ax                = axArrQdShHist(itAB);
%             bar( ax, phiCntsSh(bstSh,:) );
%             title( ax, {round(maxCrnRat,2), sum(hHist.BinCounts)} );
%         end

        % Plot histogram
        if contains(prms.measuresToPlot{itMs},'phiMax')
            if ~prms.phiHistQuadAvs
                dataForGr = mod(dataForGr+pi/4,pi*2);
                histBins  = linspace(0,2*pi,61);
                crnrLines = (0:(pi/2):(3*pi/2)) + pi/4;
            else
                dataForGr = mod(dataForGr+pi/4,pi/2);
                histBins  = linspace(0, pi/2, 16);
                crnrLines = pi/4;
            end
        end
        ax    = axArrHst( itMs, itAB );
        hHist = histogram( ax, dataForGr, histBins, 'DisplayStyle','bar','edgecolor','none' );

        % Histogram formatting
        if contains(prms.measuresToPlot{itMs},'phiMax')
            hold(ax,'on');
            for ii=crnrLines;  plot(ax, [1 1].*ii, ax.YLim, 'k:');   end
            if ~prms.phiHistQuadAvs
                tickPos = (0:(pi/2):(3*pi/2))+pi/4;
                tickLbl = {['0' char(176)],['90' char(176)],['180' char(176)],['270' char(176)]};
            else
                tickPos = [0 pi/4 pi/2];
                tickLbl = {['-45' char(176)],['0' char(176)],['45' char(176)]};
            end
            set(ax,'xTick', tickPos, 'xticklabel', tickLbl, 'xticklabelrotation',45);
        end
        ylabel(ax,'Cell count');
        xlabel(ax, prms.xLabels{itMs}, 'Interpreter','tex');
        
        % Scatter plot of phi vs d
        % - In the end I needed to do this plot in python, so all we do here is package data for export.
        dForPlot      = 17-dMaxForGr;
        dAndPhiForPy{1,itAB} = dForPlot;
        dAndPhiForPy{2,itAB} = dataForGr;

       
    end
end
% gra_multilabel( FigHst.hFig, 'col', prms.ageLabels );

save( 'dAndPhiForPy_circ.mat', 'dAndPhiForPy', '-v7');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 - analysis of ensemble shifts between BSL and CCE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2a. Define which cells, which trials to use for tunings etc.
% Get baseline scores - can be complex as can use tr1, but only if tr2 is non-BVC.
useTr1Ind        = Res.isBVC & ~Res.isBVCByTrial(:,2);   %  false(height(Res),1); % 
phiSq            = Res.phiMax(:,2);
phiSq(useTr1Ind) = Res.phiMax(useTr1Ind,1);
dSq              = Res.dMax(:,2);
dSq(useTr1Ind)   = Res.dMax(useTr1Ind,1);

% Define index of cells to use in analysis.
[uniqueDataSets, ~, dsID] = unique(Res.dataset, 'stable');
if strcmp( prms.BVCDefTrial, 'bsl' )
    valClInd = Res.isBVC & all(Res.dataInd(:,[2 5]),2);
elseif strcmp( prms.BVCDefTrial, 'cce' )
    valClInd = Res.isBVC_CCE & all(Res.dataInd(:,[2 5]),2);
elseif strcmp( prms.BVCDefTrial, 'both' )
    valClInd = Res.isBVC & Res.isBVC_CCE & all(Res.dataInd(:,[2 5]),2);
end
if ~isempty( prms.exclLongDBVCsThr )
    valClInd = valClInd & dSq <= prms.exclLongDBVCsThr;
end
if prms.exclBadCCEBVCs
    ageInBins  = ones(height(Res),1);
    for itAB=2:size(prms.ageBins,1)
        ageInBins(  Res.age>=prms.ageBins(itAB,1) & Res.age<=prms.ageBins(itAB,2) ) = itAB;
    end
    thrByAge   = prms.exclBadCCEBVCThrs( ageInBins );
    valClInd   = valClInd & Res.BVCRespMax(:,5) > thrByAge';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2b. Loop through all ensembles, get phi tunings and sq->ci shifts for later
%     analysis. Also more detailed phi & shift plots for nominated examplar
%     ensemble(s).
if ~isempty(prms.egEnsembles)
    polAxInd           = false( length(prms.egEnsembles), 3 );
    polAxInd( :, 1:2 ) = true;
    hFigEgs            = gra_multiplot( length(prms.egEnsembles), 3, 'polAx', polAxInd, 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25 );   
    axArrEgs           = getappdata(hFigEgs,'axesHandles');
    egRowCnt           = 1;
end
[cnrdPhiDiffAll, alndPhiAll, allPhiSq] = deal( cell(1,size( prms.ageBinForEns, 1 )) );
ensRVSqCi = nan(length(uniqueDataSets),3);
for itDS = 1:length(uniqueDataSets)
    dsInd   = dsID==itDS;
    dsClInd = dsInd & valClInd;
    if sum(dsClInd) < prms.minNCells
        continue;    
    end
    age       = Res.age( find(dsInd,1) );
    ageBinEns = find(  age>=prms.ageBinForEns(:,1) & age<=prms.ageBinForEns(:,2)  );
 
    phiSqCi = [phiSq(dsClInd), Res.phiMax(dsClInd,5)];
    phiDiff = circ_dist( phiSqCi(:,1), phiSqCi(:,2) );

    % Get collated across-ensemble data: centre all phi shifts around mean shift, 
    % and align all absolute circle phis by their mean shift.
    meanPhiDiff = circ_mean( phiDiff );
    cnrdPhiDiff = circ_dist( meanPhiDiff, phiDiff );
    if prms.corrEnsRots
        alndPhi   = mod( phiSqCi(:,2)+meanPhiDiff, 2*pi );
    else
        alndPhi   = phiSqCi(:,2);
    end
    cnrdPhiDiffAll{ageBinEns} = [cnrdPhiDiffAll{ageBinEns}; cnrdPhiDiff];
    alndPhiAll{ageBinEns}     = [alndPhiAll{ageBinEns}; alndPhi];
    allPhiSq{ageBinEns}       = [allPhiSq{ageBinEns}; phiSq(dsClInd)];

    ensRVSqCi(itDS,1)         = circ_r( mod(phiSq(dsClInd).*4,2*pi) );
    ensRVSqCi(itDS,2)         = circ_r( mod(Res.phiMax(dsClInd,5).*4,2*pi) );
    ensRVSqCi(itDS,3)         = ageBinEns;
    
    % Plot tuning and tuning rotation histograms for example ensemble(s). (And rate maps?).
    if any(  strcmp(prms.egEnsembles, Res.dataset{ find(dsInd,1) }  )  )
        
        % Polar histograms of phi in Sq and Ci
        T = {   phiSq(dsClInd),   Res.phiMax(dsClInd,5)  };
        for itEn=1:2
            ax        = axArrEgs(egRowCnt, itEn);
            polarhistogram(ax, T{itEn}-0.01, 60 );
            thetaticks(ax, [0 90 180 270]);
            rLab      = cell(1,length(ax.RTick));
            rLab{end} = max( ax.RTick );
            rticklabels( ax, rLab );
        
            title(ax,  ['Quad. RV = ' num2str( circ_r( mod(T{itEn}.*4,2*pi) ), '%3.2f' ) ] , 'FontWeight', 'normal' );
        end
        
        % Phi shifts Sq->Ci
        ax            = axArrEgs(egRowCnt, 3);
        histogram(ax, phiDiff, linspace(-pi,pi,60),  'DisplayStyle','bar','edgecolor','none' );
        ax.XTick      = [-pi -pi/2 0 pi/2 pi];
        ax.XTickLabel = round( ax.XTick./pi.*180 );
        hold(ax,'on');
        hLn           = plot( ax, [1 1].*meanPhiDiff, ax.YLim, 'k:' );
        hLn.LineWidth = 1.5;
        ax.YTick      = 0:max(ax.YLim);
        xlabel( ax, [char(981) ' shift Sq->Ci (' char(176) ')'] );
        ylabel( ax, 'Cell count' );
        
        egRowCnt = egRowCnt + 1;
        
        if 0
            Res.isBVCByTrial = [Res.isBVCByTrial, false( height(Res), 2), Res.isBVC_CCE];
            plotAllCellsFromTable_v2(Res,[],M,'plotFixedInd',find(dsClInd),'nRowsPerPlot',min([10 sum(dsClInd)]),'trToPlot',[1 2 5]);
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2c. Plot and analyse changes in indivdual ensemble Quad RVs across Sq and Ci.
ensRVSqCi        = ensRVSqCi( ~isnan(ensRVSqCi(:,1)), : );
adInd            = ensRVSqCi(:,3) == size( prms.ageBinForEns, 1 );
ensRVSqCi        = ensRVSqCi(:,1:2);
M                = median( ensRVSqCi, 1 );
E                = M - cat(3, prctile( ensRVSqCi, 25 ),  prctile( ensRVSqCi, 75 ) );

FigRQ   = gra_multiplot( 1, 1, 'plotsize', [4 4], 'axesborder', [1 1 1 1].*1.25 , 'structOut', 1 );  
ax      = getappdata(FigRQ.hFig,'axesHandles');
if 0
    hBar    = errorbar( ax, [1 2], M, E(:,:,1), E(:,:,2), 'bs-'  );
    set( hBar, 'linewidth', 1.5, 'markerfacecolor', 'auto' );
else
    boxplot(ax, ensRVSqCi, 'notch', 'on');
end

set( ax, 'xlim', [0.5 2.5], 'xtick', [1 2], 'xticklabel', {'Sq' 'Ci'} );
xlabel( ax, 'Environment' );
ylabel( ax, ['Ensemble quad. ' char(981) ' RV (a.u.)'] );

hold(ax,'on');
hLn                    = plot(ensRVSqCi');
[hLn(:).LineWidth]     = deal(0.5);
[hLn(adInd).Color]     = deal( [1 1 1].*0.85 );
[hLn(~adInd).Color]    = deal( [1 1 1].*0.65 );

[~,pval,~,stats] = ttest( ensRVSqCi(:,1), ensRVSqCi(:,2));
pval_wt          = signrank( ensRVSqCi(:,1), ensRVSqCi(:,2));
fprintf(1, '\nEnsemble RV Sq->Ci difference: T-test, T(%d)=%3.2f, p=%4.3f; Wilcoxon, p=%4.3f', stats.df, stats.tstat, pval, pval_wt  );
statsOut.ensRQ_p = pval_wt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2d. Plot ensemble-centered shifts, and 'de-rotated' phi, by age bin.
hFig    = gra_multiplot( 3, size( prms.ageBinForEns, 1 ), 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.4  );  
axArray = getappdata(hFig,'axesHandles');
for itAB= 1:size( prms.ageBinForEns, 1 )
    % Difference between indivdual cell shifts and mean shift (i.e. what is 
    % the within ensemble 'jitter' in shifts?).
    ax = axArray(1,itAB);
    d  = cnrdPhiDiffAll{itAB};
    histogram(ax,d,linspace(-pi,pi,60),'DisplayStyle','bar','edgecolor','none');
    hold(ax,'on');
%     plot(ax,[1 1].*-pi/4,ax.YLim,'k:');
%     plot(ax,[1 1].*pi/4,ax.YLim,'k:');
    plot(ax,[1 1].*0,ax.YLim,'k:');
    tickLabs = {'-135','-90','-45','0','45','90','135'};
    tickLabs = cellfun(@(x) [x char(176)], tickLabs, 'UniformOutput',false);
    set(ax,'xtick',[-3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4], 'xticklabel', tickLabs);
    xtickangle(ax,45);
    xlabel( ax, {[char(981) ' rotation'], ['relative to ensemble mean (' char(176) ')']});
    ylabel( ax, 'Cell Count');
    fprintf(1,'\nMedian absolute shift relative to mean, AB %d = %2.0f', itAB, circ_rad2ang(median(abs(d))));

    % Quantification of deviation of single cell phi shifts from ensemble mean phi shift.
    per45 = sum( abs(d)<=pi/4 ) / sum(~isnan(d)) * 100;
    k     = circ_kappa(d);
    fprintf(1, '\nPhi deviation from mean, Age Bin %d, %2.0f within 45deg; k=%3.2f', itAB, per45, k );

    % All phis, aligned by ensemble mean shift (supposed to be detecting 4-fold symmetry hidden by different ensemble shifts).
    ax = axArray(3,itAB);
    histogram(ax,alndPhiAll{itAB},linspace(0,2*pi,60), 'DisplayStyle','bar','edgecolor','none' );
    hold(ax,'on');
    for ii=0:(pi/2):(3*pi/2);  plot(ax,[pi/4 pi/4]+ii,ax.YLim,'k:');   end
    set(ax,'xTick',(0:(pi/2):(3*pi/2))+pi/4, 'xticklabel',{['0' char(176)],['90' char(176)],['180' char(176)],['270' char(176)]}, ...
        'xticklabelrotation',45);
    ylabel(ax,'Cell count');
    xlabel(ax, prms.xLabels{itMs}, 'Interpreter','tex');
    % Quad RV test for above.
    dQuad = mod( alndPhiAll{itAB}.*4, pi*2 );
    [p,z] = circ_rtest(dQuad);
    fprintf(1,'\nR-test (quad, post-rot) Age bin %d, z=%3.1f, p=%4.3f', itAB, z, p);
end
gra_multilabel( hFig, 'col', prms.ageLabelsEns );

% k-test for equality of K's between centred Phi Diffs for different age bins.
% i.e. is the 'jitter' of phi rotations around mean rotation different at different ages.
[p, f] = circ_ktest(cnrdPhiDiffAll{1}, cnrdPhiDiffAll{end});
fprintf(1,'\n K-Test for phi rotation diff-to-mean, AB 1 vs Adult: %4.3f, %4.3f \n', f, p);

 
