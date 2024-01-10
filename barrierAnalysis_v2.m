function [FigMS,FigSZ,statsOut] = barrierAnalysis_v2( Res, area, varargin )
% Run analyses specific to how BVCs/Border cells respond to barrier insertion.

% ANALYSIS PARAMETERS:
prms.ageBins            = [16 18; 19 21; 22 25; 26 36; 40 40];
prms.bslTrial           = 2;  % Which baseline trial to use classify the tuning of the BVC/bordC?
prms.useOnlyIsBVCTrials = 1; % if true prms.barrAnalysis.bslTrial has to also pass BVC select criterion

% What to plot %
prms.measuresToPlot     = {'barrScoreNorm' };   %
prms.groupsToPlot       = [1 2]; % Which groups (based on BVC, BordC classifiaction). For key see block 1 below.
prms.filtLongDThr       = [];  % EXCLUDE BVCs with a distance tuning longer than this.

prms.ageLabels          = {''};



% prms.doShuffle          = 0;
% prms.nShufflesPerCell   = 2000;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Basic preparatory processing %

% % load model maps
% if ~isempty( prms.modMaps)
%     Res = addModelMapsToRes( Res, prms.modMaps );
% else
%     Res = addModelMapsToRes( Res, prms.path2ModMaps );
% end

% for Moser data we want to use the first baseline trial as the second is
% a post probe trial - CHECK THAT'S CORRECT!!!!!!!!!!
if area == 2
    prms.bslTrial = 1;
    mapSz = 28;
else
    mapSz = 25;
end

% Exclude long-range tuning BVCs, if requested %
if ~isempty(prms.filtLongDThr)
    Res = Res( Res.dMax(:, prms.bslTrial ) < prms.filtLongDThr, : );
end
% Exclude high sigma zero tuning BVCs, if requested %
if ~isempty(prms.filtHighSZ)
    ind = Res.SZMaxInd(:, prms.bslTrial ) >= prms.filtHighSZ(1) & Res.SZMaxInd(:, prms.bslTrial ) <= prms.filtHighSZ(2);
    Res = Res( ind, : );
end

% exclude badly sampled data as well as data from different brain region.
Res = Res(Res.dataInd(:,prms.bslTrial) & Res.area == area, : );

if prms.useOnlyIsBVCTrials
    Res.isBVC(~Res.isBVCByTrial(:,prms.bslTrial)) = false;
    Res.isBC(~Res.isBCByTrial(:,prms.bslTrial)) = false;
end

% Combine these into group classifiers with numerical labels, so for plotting it is easy to loop over those groups of interest:
if prms.useOnlyIsBVCTrials
    Res.cellTypeGroup_1 = Res.isBVCByTrial(:,prms.bslTrial);
    Res.cellTypeGroup_2 = Res.isBCByTrial(:,prms.bslTrial);
else
    Res.cellTypeGroup_1 = Res.isBVC;
    Res.cellTypeGroup_2 = Res.isBC;
end
Res.cellTypeGroup_3 = Res.cellTypeGroup_1 & Res.cellTypeGroup_2;
Res.cellTypeGroup_4 = Res.cellTypeGroup_1 & ~Res.cellTypeGroup_2;
Res.cellTypeGroup_5 = ~Res.cellTypeGroup_1 & Res.cellTypeGroup_2;
Res.cellTypeGroup_6 = ~Res.isBVC;  % This is deliberate - want this group to always be cells that are never a BVC across baselines.


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) For each cell, get the dir tuning: get the barrier orientation orthogonal to the prefered tuning and also work out which side of the barrier the 'doubled' field should fall.
   
% 2a. Need the preferred wall as defined by BVC analysis.
% Note new version of this - previously including only phi tunings within a 45 deg window of cardinal dirs.
% Now only exlcudes most extreme 'ambiguous' cases (within one ang bin of tuning towards box corner).
phiSh45           = mod( Res.phiMax( :, prms.bslTrial ) + (pi/4), 2*pi );  % Shift for convenience (quadrants start at 0 not 7pi/4).
quads             = [0 pi/2 pi 3*pi/2 pi*2];
BVCQuad           = nan(size(phiSh45));
crnExcAng         = 18  /180*pi;   % Exclude BVCs with corner tunings (as ambiguous which barrier to use). This variable sets half-width of exclusion window around 45 deg.
for itQd = 1:4
    BVCQuad(   phiSh45>=(quads(itQd)+crnExcAng) & phiSh45<(quads(itQd+1)-crnExcAng)   ) = itQd;
end

% Re-order quad flags so that the match the 'borderScoreWall' convention: 1=N, 2=W, 3=S, 4=E - check that N-S is correct as based on image coords (old convention: 1=W, 2=E, 3=N, 4=S)
reOrder                  = [4 3 2 1];            % Matlab won't let me declare and index on the same line? [4 3 2 1]; [2 4 1 3] - for older data before new borderScore calc is introduced (pre Jan-15 2019)  
BVCQuad(~isnan(BVCQuad)) = reOrder( BVCQuad(~isnan(BVCQuad)) )';
Res.BVCQuad              = BVCQuad ; 

% 2b. On the basis of this, create a logical index to select trials.
[Res.orthBarrTrBVC, Res.orthBarrTrBoC] = deal( false( size(Res,1), 2 ) );
orthBarrInd = false(height(Res),1);
orthBarrInd(~isnan(Res.BVCQuad)) = ~mod(Res.BVCQuad(~isnan(Res.BVCQuad)),2);
Res.orthBarrTrBVC( orthBarrInd & ~isnan(Res.BVCQuad), 1 )  = true;   % First of barrier trials in standard sequence is N-S, so orthog to cells with E or W wall tunings.
Res.orthBarrTrBVC( ~orthBarrInd & ~isnan(Res.BVCQuad), 2 ) = true;   % Second of barrier trials in standard sequence is W-E, so orthog to cells with N or S wall tunings.

bordScoreWall                           = Res.borderScoreWall( :, prms.bslTrial );
orthBarrIndBoC                          = false(size(bordScoreWall));   % To preserve cell indexing order, we also need to add an entry for when the border score is NaN (obviously which actual trial is irrelevant).
orthBarrIndBoC(~isnan(bordScoreWall))   = ~mod(bordScoreWall(~isnan(bordScoreWall)),2);
Res.orthBarrTrBoC( orthBarrIndBoC, 1 )  = true;   % The same as the two rows above, this time for border cell 'mainWall'.
Res.orthBarrTrBoC( ~orthBarrIndBoC, 2 ) = true;   %   ..


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (4) First plot Moser style barrier scores
FigMS  = gra_multiplot( 1, length(prms.measuresToPlot), 'plotsize', [ 2 2 ], 'axesborder', [1 1 1 1].*1.25, 'structOut', 1 );   
axArrMS = FigMS.axArr;
FigSZ  = gra_multiplot( 1, length(prms.measuresToPlot), 'plotsize', [2 2], 'axesborder', [1 1 1 1].*1.25,  'structOut', 1 );   
axArrSZ = FigSZ.axArr;

cellIDs = cell( 1, size(prms.ageBins,1) );  % For export of cellIDs for representative rate maps.
for itMs = 1:length(prms.measuresToPlot)

    % Preallocate mean and error for group/age bins.
    [M,E,M2,E2,M3,E3] = deal(nan(length(prms.groupsToPlot),size(prms.ageBins,1)));
    d4Anova           = []; 
    gr4Anova          = cell(3,1);
    for itGr = 1:length(prms.groupsToPlot)
        
        % Get a) the cells to plot (aligned EW or NS) - this is 'validCells'
        %     b) the correct barrier trial for the cell, depending on orientation - 'trSelInd'
        %     c) the orientation of the *barrier side* where new field is expected - 'prefWallInd' - values 1,2,3,4.
        if any( prms.groupsToPlot(itGr)==[1 3 4 6] )
            % This is for BVCs
            validCells  = ~isnan(Res.BVCQuad);
            trSelInd    = Res.orthBarrTrBVC(validCells,:);
            prefWallInd = Res.BVCQuad(validCells);
        elseif any( prms.groupsToPlot(itGr)==[2 5] )
            % This is Border cells - different scores to classify relevant wall
            validCells  = ~isnan( Res.borderScoreWall(:,prms.bslTrial) );
            trSelInd    = Res.orthBarrTrBoC(validCells,:);
            prefWallInd = Res.borderScoreWall(validCells, prms.bslTrial );
        end

        % Work out which barrier trial, and which side of barrier, should be doubling.
        % Following this get the actual rate scores on the expected high and expected low sides - 'dataHigh' and 'dataLow'
        zoneInd           = false( size(prefWallInd,1), 2 );
        zoneInd( prefWallInd==4 | prefWallInd==3, 1 ) = true;  % The first entry in the vector is for the 'lower' zone (lower in matlab indicing). Congruent with fields on E and S walls, so when wall inds are 4 or 1.
        zoneInd( prefWallInd==1 | prefWallInd==2, 2 ) = true;  % The 2nd entry in the vector is for the 'upper' zone (greater in matlab indicing). Congruent with fields on W and N walls, so when wall inds are 3 or 2.
        dataBothTrsBothZones = Res.( prms.measuresToPlot{itMs} )(validCells,3:4)'; 
        dataOrthTrBothZones  = cell2mat( dataBothTrsBothZones( trSelInd' ) )';
        dataHigh             = dataOrthTrBothZones( zoneInd' );  %  
        dataLow              = dataOrthTrBothZones( ~zoneInd' );

        % Get barrier scores for barrier parallel to tuning, as a control for responses to orthogonal.
        paraBarrTrInd        = ~trSelInd;
        dataParaTrBothZones  = cell2mat( dataBothTrsBothZones( paraBarrTrInd' ) )';

        % Debug
        if strcmp(prms.measuresToPlot{itMs}, 'barrScoreNorm') 
            debugData    = dataHigh; 
            debugCellInd = validCells;
        end
        
        % For each age bin, get raw data and plot histogram. Collect mean and error for later plotting.
        for itAB=1:size(prms.ageBins,1)  % itAB=iterator for age group
            ageInd   = Res.age(validCells) >= prms.ageBins(itAB,1) & Res.age(validCells) <=prms.ageBins(itAB,2);
            groupInd = Res.( ['cellTypeGroup_' num2str( prms.groupsToPlot(itGr) )] )(validCells);
            SZValCls = Res.SZMaxInd(validCells,prms.bslTrial);
            SZForGr  = SZValCls( ageInd & groupInd );
            if sum(ageInd & groupInd)==0;  continue;  end

            % High firing rate zone
            d               = dataHigh( ageInd & groupInd );
            M(itGr,itAB)    = nanmean(d);
            E(itGr,itAB)    = nanstd(d) / sqrt(sum(~isnan(d)));
            % Low firing rate zone.
            d2               = dataLow( ageInd & groupInd );
            M2(itGr,itAB)    = nanmean(d2);
            E2(itGr,itAB)    = nanstd(d2) / sqrt(sum(~isnan(d2))); 
            % Parallel barrier trial (both zones).
            d3               = dataParaTrBothZones( :, ageInd & groupInd );
            d3               = d3( ~isnan(d3) );
            randSampInd      = randi(length(d3),[sum(~isnan(d)), 1000]);
            randSampMeans    = mean( d3(randSampInd), 1 );
            M3(itGr,itAB)    = prctile(randSampMeans,95);
            E3(itGr,itAB)    = nan; % 1.64 * nanstd(d3) / sqrt(sum(~isnan(d3)));  % nanstd(d3) / sqrt(sum(~isnan(d3))); 
            
            % How many cells actually increase firing on barrier distal side?
            if itGr==1
                fprintf( '\n Perc Distal>0: %3.2f', sum(d>0)*100/length(d) );
            end

            % Get cellIDs for a closer look at the data:
            if prms.groupsToPlot(itGr)==1 && itMs==1 && itAB==3
                cellIDs_gr       = Res.cellID(validCells);
                cellIDs_gr       = cellIDs_gr( ageInd & groupInd );
                repScInd         = [1 1].*M(itGr,itAB) + [-3 3].*E(itGr,itAB);
%                 scMtchInd        = d>=repScInd(1) & d<=repScInd(2); 
                scMtchInd        = SZForGr>=3;
                cellIDs{itAB}    = [ cellIDs_gr( scMtchInd ) num2cell( d(scMtchInd) ) ];
            end
            
            % stats
            d4Anova     = [d4Anova; [d d2]];
            gr4Anova{1} = [gr4Anova{1}; itAB.*ones(length(d),1)];
            gr4Anova{2} = [gr4Anova{2}; itGr.*ones(length(d),1)];
            gr4Anova{3} = [gr4Anova{3}; SZForGr];

        end
    end  


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final Manuscript plots.
    ax               = axArrMS(itMs);  
    if 0
        axes(ax);
        % For MS plots we want barrier high and low side, only for BVCs.
        msM              = [M(1,:); M2(1,:)];
        msE              = [E(1,:); E2(1,:)];
        hBr              = barwitherr(msE',msM');
        hBr(2).FaceColor = [0.6 0.6 1];
        hBr(1).FaceColor = [0 0.1 1];
        set(ax,'xlim',[0 size(prms.ageBins,1)+0.5],'xticklabels',prms.ageLabels,'xtick',1:size(prms.ageBins,1),'XTickLabelRotation',45,'ylim',prms.yLim{itMs},'box','off' );
    else
        BVCInd           = gr4Anova{2}==1;
        d4Box            = [d4Anova(BVCInd,1); d4Anova(BVCInd,2)];
        gr4Box{2}        = [ ones(sum(BVCInd),1); ones(sum(BVCInd),1).*2];
        gr4Box{1}        = repmat( gr4Anova{1}(BVCInd), 2, 1);
        hBP              = boxplot(ax,  d4Box, gr4Box, 'Notch','on','Colors', [0 0.1 1; 0.6 0.6 1], 'Symbol', '',...
                               'FactorGap', [10 0], 'boxstyle', 'outline', 'width',0.7,'LabelVerbosity','majorminor',...
                               'labels', {'','P16-18','','P19-21','','P22-25','','Adult'} );
        set(hBP, 'LineWidth', 1);
        hold( ax, 'on' );
        plot(ax, ax.XLim, [ 0 0 ], 'k-' );
    end    
    xlabel(ax,'Age');
    ylabel(ax, prms.yLabelsMS{itMs})

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot results split by (baseline) sigma zero. Works from 'data4Anova' data.
    for itAB = 4 % 1:size(prms.ageBins,1)
        % One plot for each age bin.
        M = nan(2,4);    E = nan(2,4);
        for itSZ = 1:4
            for itSd = 1:2
                d            = d4Anova( gr4Anova{1}==itAB & gr4Anova{3}==itSZ, itSd );
                M(itSd,itSZ) = nanmean(d);
                E(itSd,itSZ) = nanstd(d) / sqrt(sum(~isnan(d)));
            end
        end
        
        % Plot

        ax = axArrSZ(1,itMs);  
        if 0
            axes(ax);
            hBr2              = barwitherr(E',M');
            hBr2(2).FaceColor = [0.6 0.6 1];
            hBr2(1).FaceColor = [0 0.1 1];
            set(ax,'xticklabels',Res.Properties.UserData.szList,'XTickLabelRotation',45,'ylim',prms.yLim{itMs},'box','off' );
        else
            BVCInd           = gr4Anova{2}==1 & gr4Anova{1}==itAB;
            d4Box            = [d4Anova(BVCInd,1); d4Anova(BVCInd,2)];
            gr4Box{2}        = [ ones(sum(BVCInd),1); ones(sum(BVCInd),1).*2];
            gr4Box{1}        = repmat( gr4Anova{3}(BVCInd), 2, 1);
            hBP              = boxplot(ax,  d4Box, gr4Box, 'Notch','on','Colors', [0 0.1 1; 0.6 0.6 1], 'Symbol', '',...
                                   'FactorGap', [10 0], 'boxstyle', 'outline', 'width',0.7,'LabelVerbosity','majorminor',...
                                   'labels', {'','6.2','','12.2','','20.2','','30.2'} );
            set(hBP, 'LineWidth', 0.8);
            hold( ax, 'on' );
            plot(ax, ax.XLim, [ 0 0 ], 'k-' );
        end
        xlabel(ax,[char(963) '_0 (cm)'], 'interpreter', 'tex');
        ylabel(ax, prms.yLabelsMS{itMs});
        gra_multilabel(FigSZ.hFig,'row',prms.ageLabels);



        % Stats - ANOVA for prox firing, factor = sig zero
        anovaInd      = gr4Anova{2}==1 & gr4Anova{1}==itAB;  % Select just BVCs and one age bin at a time.
        [p,tbl,stats] = anova1(d4Anova(anovaInd,2),gr4Anova{3}(anovaInd),'off');
        fprintf('\nSZ (AB=%i): F(%i,%i)=%3.2f; p=%6.5f', itAB, tbl{2,3},tbl{end-1,3},tbl{2,5},p(1));
        statsOut.SZSplit_ANOVA{itAB} = [ tbl{2,3},tbl{end-1,3},tbl{2,5},p(1) ];
        % Stats - t-test all prox against 0
        fprintf( 1, '. \t T-test prox vs 0, AB=%i', itAB );
        for itSZ=1:4
            [~,p,~,S] = ttest( d4Anova( gr4Anova{1}==itAB & gr4Anova{3} == itSZ, 2 ), 0, 'tail', 'left' );
            fprintf( '\t SZ=%i, t(%d)=%2.3f, p=%4.3f', itSZ, S.df, S.tstat, p);
            statsOut.ttest_0_SZ{itAB,itSZ} = [S.df, S.tstat, p];
        end
    end
    fprintf(1,'\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Export data for SPSS Age*Side ANOVA (so can use interaction+SME tests)
    % Export this for BVCs only - assuming these are group 1.
    BVCInd = gr4Anova{2}==1;
    d4SPSS = [ d4Anova( BVCInd, : ), gr4Anova{1}(BVCInd) ];
    assignin('base',['d4SPSS_' prms.measuresToPlot{itMs}],d4SPSS);
    % Age*barr side ANOVA, for all ages.
    tbl = simple_mixed_anova( d4Anova( BVCInd, : ),  gr4Anova{1}(BVCInd), {'side'}, {'age'} );
    statsOut.ANOVA = [tbl.DF(5), tbl.DF(6), tbl.F(5), tbl.pValue(5)];
    % Test one cell type against another, for one barrier side at a time.
    for itBS=1:2
        [p,tbl,stats] = anovan(d4Anova(:,itBS),gr4Anova(1:2),'display','off','model','full','varnames',{'AB','isBVC'});
        fprintf('\n\n %s, Barrside %d', prms.measuresToPlot{itMs}, itBS)
        fprintf('Age: F(%i,%i)=%.2f; p=%.5f\n',tbl{2,3},tbl{end-1,3},tbl{2,6},p(1));
        fprintf('Factor Group: F(%i,%i)=%.2f; p=%.5f\n',tbl{3,3},tbl{end-1,3},tbl{3,6},p(2));
        fprintf('Factor Age*Group: F(%i,%i)=%.2f; p=%.5f\n',tbl{4,3},tbl{end-1,3},tbl{4,6},p(3));
        [PWComp,m,h,gnames] = multcompare(stats,'Dimension',[1 2],'Display','off');
        prs2Show = [1 5; 2 6; 3 7; 4 8];
        for itPr = 1:size(prs2Show,1)
            if 1
                % Option 1 for pairwise is t-test, and adjust with sidak (p thr for 4*tests = 0.0127)
                ind1   = gr4Anova{1}==itPr & gr4Anova{2}==1;
                ind2   = gr4Anova{1}==itPr & gr4Anova{2}==2;
                [~,p]  = ttest2( d4Anova(ind1,itBS), d4Anova(ind2,itBS) );
            else
                % Option 2 is pairwise HSD by multcompare. Beware, p's are counter-intuitively low, I think
                % because inclusion of adult data makes pooled variance very high. T-tests more powerful, even if
                % theoretically more conservative.
                ind = PWComp(:,1)==prs2Show(itPr,1) & PWComp(:,2)==prs2Show(itPr,2);
                p   = PWComp(ind,6);
            end
            fprintf( '\tAB%d, p=%4.3f', itPr, p);
        end
    end
    % Test proximal rates versus zero.
    for itAB=1:size(prms.ageBins,1)
        [~,p,~,S] = ttest( d4Anova( gr4Anova{1}==itAB, 2 ), 0, 'tail', 'left' );
        fprintf( '\n T-test prox vs 0, AB%d, t(%d)=%2.3f, p=%4.3f', itAB, S.df, S.tstat, p);
        statsOut.ttest_0{itAB} = [S.df, S.tstat, p];
    end

end

% Export to base workspace cell IDs of cells with representative barrier scores
assignin('base',"cellIDs",cellIDs);
