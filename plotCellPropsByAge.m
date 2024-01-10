function [Fig, statsOut] = plotCellPropsByAge( Res, area, varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%% Params
prms.cellType  = 'BVC';
prms.ageBins   = [16 18; 19 21; 22 25; 26 36; 40 40];  %  default: [16 18; 19 21; 22 25; 26 36; 40 40];     [15 18; 19 21; 22 24; 25 27; 28 30; 31 33; 34 36; 40 40];  %   [15 17; 18 19; 20 21; 22 25; 26 36; 40 40];   %    
prms.ageLabels = {''};
% Plotting %
prms.bslTrToPlot          = 1:2;
prms.plotOnlyIsBVCTrials  = 1;
prms.splitVar             = []; % 'dMax';   %   'env'    %   'cellType' 'cellXarea'  % How to split data for >1 line per plot. 'env'=split by bsl/barr/CCE, 'cellType'=if plotting BVCs, those cells that are also BC, or vice versa.
prms.splitThr             = [];
prms.nEnvToPlot           = 1;  % 1=baseline only, 2=bsl + barrier, 3=basl+barrier+CCE. Only effective if prms.splitVar='env', otherwise always just plot baselines.
prms.measuresToPlot       = { 'SIAdSm', 'intraStabAdSm', 'interStabAdSm','borderScoreAdSm' }; % ,'intraStabAdSm', 'interStabAdSm', , 'interStabHT', 'BVCRespMax', 'dMax' ,'diffToBestD', 'diffToBestPhi','rAtBestBSLPhiD'};
prms.yLimMax              = {0.8,       0.9,              0.9,           0.8,           1,             1,            0.8,             0.8      , [],            [],             []};

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

%%
% Plot spatial props %

Fig = gra_multiplot( length(prms.measuresToPlot), 1, 'plotsize', [3 2],'axesborder', [1 1 1 1].*1.5, 'structOut', 1 );   axArray=Fig.axArr;

if ~isempty( prms.splitVar )
    nLineGroups = 2;
else
    nLineGroups = prms.nEnvToPlot;
end

for ii=1:length(prms.measuresToPlot)   
    
    % Apply filters for poor data sampling. For inter-trial measures (which are inter-baseline), make this such that both baselines must be well-sampled. %
    if strncmp( prms.measuresToPlot{ii}, 'inter', 5 )
        sampleFilterForBsl = all( Res.dataInd(:,1:2), 2 );
        Res.( prms.measuresToPlot{ii} )( ~sampleFilterForBsl ) = nan;
    elseif strcmp(prms.measuresToPlot{ii},'distPhi')
        Res.phiMax( ~Res.dataInd ) = nan;
    elseif strcmp(prms.measuresToPlot{ii},'distD')  
        Res.dMax( ~Res.dataInd ) = nan;
    else
        Res.( prms.measuresToPlot{ii} )( ~Res.dataInd ) = nan;
    end
    
    % Get the Means and SEMs within each group %
    [M,E] = deal( nan(nLineGroups,size(prms.ageBins,1)) );
    [dataForStats, groupIndForStats] = deal([]);
    alpha = []; idx = []; deltaD = []; groupD = []; 
    for jj=1:size(prms.ageBins,1)
        
      
        switch prms.cellType
            case 'BVC'
                if strcmp(prms.splitVar, 'cellType' )
                    indForDatapoint = (Res.isBVC | Res.isBC) & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2) & Res.area == area;
                    isCellTypeByTrial = Res.isBVCByTrial(indForDatapoint,:) | Res.isBCByTrial(indForDatapoint,:);
                elseif strcmp(prms.splitVar, 'area' )
                    indForDatapoint = Res.isBVC & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2);
                    isCellTypeByTrial = Res.isBVCByTrial(indForDatapoint,:);
                else
                    indForDatapoint = Res.isBVC & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2) & Res.area == area;
                    isCellTypeByTrial = Res.isBVCByTrial(indForDatapoint,:);
                end
            case 'BC'
                if strcmp(prms.splitVar, 'cellType' )
                    indForDatapoint = (Res.isBVC | Res.isBC) & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2) & Res.area == area;
%                     isCellTypeByTrial = Res.isBVCByTrial(indForDatapoint,:) | Res.isBCByTrial(indForDatapoint,:);
                elseif strcmp(prms.splitVar, 'area' )
                    indForDatapoint = Res.isBC & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2);
                    isCellTypeByTrial = Res.isBCByTrial(indForDatapoint,:);
                else
                    indForDatapoint = Res.isBC & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2) & Res.area == area;
                    isCellTypeByTrial = Res.isBCByTrial(indForDatapoint,:);
                end
             case '~BVC'
                indForDatapoint = ~Res.isBVC & Res.age>=prms.ageBins(jj,1) & Res.age<=prms.ageBins(jj,2) & Res.area == area;
                isCellTypeByTrial = ~Res.isBVCByTrial(indForDatapoint,:);
        end
        
        if strcmp(prms.measuresToPlot{ii},'distPhi')
            d = circ_dist(Res.phiMax(indForDatapoint,1),Res.phiMax(indForDatapoint,2));
            alpha = [alpha; d]; idx = [idx;jj*ones(length(d),1)];
        elseif strcmp(prms.measuresToPlot{ii},'distD')
            d = abs(Res.dMax(indForDatapoint,1)-Res.dMax(indForDatapoint,2));
            deltaD = [deltaD; d]; groupD = [groupD;jj*ones(length(d),1)];
        else
            d  = Res.( prms.measuresToPlot{ii} )(indForDatapoint,:);
        end
        
        if strncmp( prms.measuresToPlot{ii}, 'inter', 5 )
            d = [d nan(size(d,1), 4)];    % Need to pad inter-trial measures (which are always nCol=1, just compare baselines), so that the indexing in the section below works.
        end
        
        % Sort the data into different means by env type %
        if prms.nEnvToPlot==1
            % If requested, only plot those trials in which the cell is actually classified as BVC/BC. If only allowed this for baseline only plot (prms.nGroupToPlot=1),
            % as otherwise plots wouldn't be equivalent between baseline and probe.
            if prms.plotOnlyIsBVCTrials
                if ~isempty(prms.splitVar)
%                     ind1               = isCellTypeByTrial &  Res.isBVCByTrial(indForDatapoint,:);
%                     ind2               = isCellTypeByTrial &  Res.isBCByTrial(indForDatapoint,:);
                    if strcmp(prms.splitVar, 'cellType' ) 
                        ind1               = Res.isBVCByTrial(indForDatapoint,:);
                        ind2               = Res.isBCByTrial(indForDatapoint,:);
                    elseif strcmp(prms.splitVar, 'area' ) 
                        ind1               = isCellTypeByTrial & Res.area(indForDatapoint) == 1;
                        ind2               = isCellTypeByTrial & Res.area(indForDatapoint) == 2;
                    end
                    [temp_d1, temp_d2] = deal(d);
                    if strncmp( prms.measuresToPlot{ii}, 'inter', 5 )
                        temp_d1( sum(ind1,2)<1, : ) = nan;
                        temp_d2( sum(ind2,2)<1, : ) = nan;
                    else
                        temp_d1( ~ind1 ) = nan;
                        temp_d2( ~ind2 ) = nan;
                    end
                    if strcmp(prms.measuresToPlot{ii}(1:2),'di')
%                         d = [circ_mean(temp_d1(:,prms.bslTrToPlot ),[],2) circ_mean(temp_d2(:,prms.bslTrToPlot ),[],2)];
                    else
                        d = [nanmean(temp_d1(:,prms.bslTrToPlot ),2) nanmean(temp_d2(:,prms.bslTrToPlot ),2)];
                    end
                else
                    if strncmp( prms.measuresToPlot{ii}, 'inter', 5 )
                        d( sum(isCellTypeByTrial,2)<1, : ) = nan;
                    elseif strcmp(prms.measuresToPlot{ii}(1:2),'di')
                    else
                        d( ~isCellTypeByTrial ) = nan;
                    end
                    
                    if strcmp(prms.measuresToPlot{ii}(1:2),'di')
%                         d = circ_mean(d(:,prms.bslTrToPlot ),[],2);
                    else
                        d = nanmean(d(:,prms.bslTrToPlot ),2);
                    end
                end
            else
                if strcmp(prms.measuresToPlot{ii}(1:2),'di')
%                     d = circ_mean(d(:,prms.bslTrToPlot ),[],2);
                else
                    d = nanmean(d(:,prms.bslTrToPlot ),2);
                end
            end
        elseif prms.nEnvToPlot==2
            if strcmp(prms.measuresToPlot{ii}(1:2),'di')
%                 d = [circ_mean(d(:,prms.bslTrToPlot ),[],2); circ_mean(d(:,3:4),[],2)];
            else
                d = [nanmean(d(:,prms.bslTrToPlot ),2), nanmean(d(:,3:4),2)];  % When there are 4Col scores, take the mean for bsl and mean for barr.   
            end
        elseif prms.nEnvToPlot==3
            if strcmp(prms.measuresToPlot{ii}(1:2),'di')
%                 d = [circ_mean(d(:,prms.bslTrToPlot ),[],2); circ_mean(d(:,3:4),[],2), d(:,5)];
            else
                d = [nanmean(d(:,prms.bslTrToPlot ),2), nanmean(d(:,3:4),2), d(:,5)];  % When there are 4Col scores, take the mean for bsl and mean for barr. 
            end
        end
        
        % Split the data by belonging to the alternate cell type, or an arbitrary variable (presumes you are only plotting bsl, i.e. prms.nEnvToPlot = 1) %
        if ~isempty( prms.splitVar )
            if strcmp(prms.splitVar, 'cellType' )  %%%%%% DOESNT WORK ATM WITH prms.plotOnlyIsBVCTrials=0 %%%%%%%%%%%%%%%%%%
                
                splitInd1 = Res.isBVC(indForDatapoint);
                splitInd2 = Res.isBC(indForDatapoint);
%                 d(~splitInd,1) = nan; %%% LM: This seems wrong to do for cell type as we can't just arbitrarily count some cells to one group but n0t the other
                d(~splitInd1,1) = nan; 
                d(~splitInd2,2) = nan;  % Instead only remove from cells that are 
            elseif strcmp(prms.splitVar, 'area' )  
                % already split!
            else
                splitInd = Res.( prms.splitVar )(indForDatapoint, prms.thrTrials) <= prms.splitThr;
                d = repmat(d, 1, 2);
                d(~splitInd,1) = nan;
                d(splitInd,2) = nan;
            end
        end
        
        % Get a record of the final raw datapoints contributing to mean and SEM, for later stats %
        dataForStats     = [dataForStats; d(:)]; %#ok<AGROW>   % d(:) is quick fix to avoid having to deal with d that is split by splitVar: really need to set up a second group var.
        groupIndForStats = [groupIndForStats; ones(numel(d),1).*jj]; %#ok<AGROW>
        % Finally, get the mean and SEM for the plot %
        for kk=1:nLineGroups
            if strcmp(prms.measuresToPlot{ii},'distPhi')
                M(kk,jj) = abs(circ_mean( d(:,kk) ) );
                E(kk,jj) = circ_std( d(:,kk) ) / sqrt( sum(~isnan( d(:,kk)  )) ); 
            else
                M(kk,jj) = nanmean( d(:,kk) );
                E(kk,jj) = nanstd( d(:,kk)  ) / sqrt( sum(~isnan( d(:,kk)  )) );   
            end
        end
    end
    
    offset = [-1 0 1].*0.1;  colours = {'b','r','m'};
    
%     [p,tbl,stats] = anova1(dataForStats,groupIndForStats);
%     c = multcompare(stats);
   
    % Plot (errorbar plot of Mean +/- SEM) %
    ax = axArray(ii);
    if 0
        for jj=1:nLineGroups
            errorbar(ax, (1:size(M,2))+offset(jj), M(jj,:), E(jj,:), [colours{jj} 'o-'],'markerfacecolor',colours{jj},'markeredgecolor','none','linewidth',1,'markersize',1);
            hold(ax,'on');
        end
        set(ax,'ylim',[0 max(get(ax,'ylim'))], 'xlim', [0.5 size(prms.ageBins,1)+0.5], 'xtick', 1:size(M,2), 'xticklabel', prms.ageLabels, 'XTickLabelRotation', 45, 'box', 'off');
    else
        % Create a fake second grouping variable, that splits ad vs pup, so that ad box is spaced a bit away from pup boxes.
        gr4Box                = groupIndForStats;
        fakeGrAd              = ones(size(gr4Box));
        fakeGrAd( gr4Box==4 ) = 2;
        gr4Box(gr4Box==4)     = 1;
        % PLot
        hBP  = boxplot(ax,  dataForStats, {fakeGrAd, gr4Box}, 'Notch','on','Colors', 'b', 'Symbol', '',...
                  'FactorGap', [10 0], 'boxstyle', 'outline', 'width',0.6,'LabelVerbosity','majorminor',...
                  'labels', prms.ageLabels );
        set(hBP, 'LineWidth', 0.8);
        if prms.yLimMax{ii}(1) < 0
            hold( ax, 'on' );
            plot(ax, ax.XLim, [ 0 0 ], 'k-' );
        end
        set(ax,'box','off');
    end
    if ~isempty(prms.yLimMax{ii});  set(ax,'ylim',prms.yLimMax{ii});    end
    ylabel(ax,prms.yLabels{ii},'Interpreter','tex');
    
    % Run ANOVA and add results to plot %
    if isempty( strfind( prms.splitVar, 'cell' ) )
        [pVal,T,STATS] = anova1( dataForStats, groupIndForStats, 'off' );
        C              = multcompare(STATS,'display','off');
        fprintf(1,'\n%s: F(%d,%d)=%3.2f, p=%4.3f; \t', prms.measuresToPlot{ii}, T{2,3}, T{3,3}, T{2,5}, T{2,6});
        fprintf(1,'\nAd vs P25: p=%4.3f, P16 vs P25: p=%4.3f', C(6,6), C(2,6));
        % Export stats for further use:
        statsOut(ii).F          = T{2,5};
        statsOut(ii).p          = T{2,6};
        statsOut(ii).df         = [T{2,3}, T{3,3}];
        statsOut(ii).p_pairwise = [C(6,6), C(2,6)];


%         if ~isempty(prms.yLimMax{ii})
%             text( 0.5, min(prms.yLimMax{ii}), ['p=' num2str(pVal,'%4.3f')], 'fontsize', 10, 'parent', ax, 'verticalalignment', 'bottom' );
%         else
%             text( 0.5, 0, ['p=' num2str(pVal,'%4.3f')], 'fontsize', 10, 'parent', ax, 'verticalalignment', 'bottom' );
%         end
    end
    
    
end
% gra_multilabel(Fig.hFig, 'col', prms.measuresToPlot);
% if prms.area==1;  titleStr='Sub';  else   titleStr='mEC';    end
% gra_multilabel( hFig, 'title', titleStr);

end

