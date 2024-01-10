function [ ResBVC ] = getCellTypeProportions( Res, area, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prms.ageBins       = [16 18; 19 21; 22 25;  40 40];  %  default: [16 18; 19 21; 22 25; 26 36; 40 40];     [15 18; 19 21; 22 24; 25 27; 28 30; 31 33; 34 36; 40 40];  %   [15 17; 18 19; 20 21; 22 25; 26 36; 40 40];   %    
prms.ageLabels     = {'P16-18','P19-21','P22-25','Adult'};
prms.axis          = [];
prms.cellType2Plot = 'isBVC';
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ----------------------------------------------------------------------------------


% Plot N BVCs (percentage) %
[N,D,P,E] = deal(  zeros(size(prms.ageBins,1), 1 )   );
for ii=1:size(prms.ageBins,1)
    
    ageInd          = Res.age>=prms.ageBins(ii,1) & Res.age<=prms.ageBins(ii,2);
    indForDatapoint =  ageInd & any(Res.dataInd(:,1:2),2);

    N(ii,1)  = sum(Res.(prms.cellType2Plot)( indForDatapoint & Res.area==area) );
%     N(ii,2)  = sum(Res.isBC( indForDatapoint & Res.area==area) );
%     N(ii,1) = sum( any( Res.isBVCByTrial(:,1:2) & Res.dataInd(:,1:2) & Res.area==area & ageInd, 2 ) );
    
    D(ii,1) = sum(indForDatapoint & Res.area==area);
    P(ii,:) = N(ii,:)*100 ./ D(ii,:);  
    
end
disp(N);

% [p, chi2stat] = chi2TestOfIndependence(N(:,1), D(:,1));
% for i = [1 2 3 4 5] 
% [ Z(i), p(i) ] = do_Ztest( N(i,1)/D(i,1), N(i,2)/D(i,2), D(i,1), D(i,2), 2 )
% end

pVal = [NaN;NaN];
for i = 1

    % Confidence interval for nSigEv (ents) %
    pHat      = P(:,i) ./ 100;   % Zar Eq 24.63 - convert % to prob.
    qHat      = 1 - pHat;
    E(:,i)    = 1.96  .*   sqrt(   ( pHat.*qHat ) ./ D(:,i)  )  .*  100   ./   2;


    % Z-test of proportions, adult versus pup (a la Bjerknes 2014) or between proportions of 'splitvar' %
    % Change test groups by changing index definition below:
    if 0
        testInd1 = prms.ageBins(:,1) == 40;
        testInd2 = prms.ageBins(:,1) >=22 & prms.ageBins(:,1) <= 25;
    else
        testInd1 = prms.ageBins(:,1) <=18;
        testInd2 = prms.ageBins(:,1) >=22 & prms.ageBins(:,1) <= 25;
    end
    nBVC  = [ sum( N(testInd2,i) ),  nansum( N(testInd1,i) ) ];
    nCell = [ sum( D(testInd2,i) ),  nansum( D(testInd1,i) ) ];
    pHat  = nBVC ./ nCell;
    pBar  = sum( nBVC, 2 )  ./ sum( nCell, 2 );   % Eq 24.50
    qBar  = 1 - pBar;
    Z     = (pHat( 1 ) - pHat( 2 ))    /     sqrt(    ((pBar*qBar)/nCell( 1 ))     +     ((pBar*qBar)/nCell( 2 ))    );    % Eq 24.49
    pVal(i)  = (1 - normcdf(abs(Z),0,1)).*2; % times two of the raw cdf is two-tailed probability
    fprintf(1,'\n Z-test: Z=%3.2f, p=%4.3f \n', Z, pVal(i));

end


% plot
if isempty(prms.axis)
    hFig         = gra_multiplot( 1, 1, 'plotsize', [2.25 2],'axesborder', [1 1 1 1].*1.5 );   axArray=getappdata(hFig,'axesHandles');
    axes(axArray(1));
else
    axes( prms.axis );
end
if 1
    [hBr, hErr]  = barwitherr(E,P); % SEM is 95% confidence intervals
    hBr.BarWidth = 0.6;
else
    errorbar( 1:length(P), P, E );
end

% Format plot %
ylim = [0 50];
set(gca,'xlim', [0 size(prms.ageBins,1)+1], 'xtick', 1:size(P,1), 'xticklabel', prms.ageLabels, 'XTickLabelRotation', 45, 'ylim', ylim,'ytick', 0:10:50,'box','off');
ylabel('% BVCs'); xlabel('Age');

% False postive lines 
% CAREFUL - this was added late in MS drafting, thr levels are hard-coded.
falPosRate = 1 - ((1-(0.01*0.75))^2);          % Assuming that 99% BVC and 75% SI thrs are independent, & every cell has independent chance of passing on each of two bsl trials.
sigLev95   = binoinv( 0.95, D, falPosRate ) ./ D .*100;   
hold on;
for itAB=1:4
    plot( [-0.4 0.4]+itAB, [1 1].*sigLev95(itAB), 'k:' )
end

% format plot
% str = strcat(cellstr(num2str(N(:,1))), ' /',{' '}, cellstr(num2str(N(:,2))) );
% fprintf('%s\n','===== Proportions of BVCs & BCs ======');
% str2 = strcat(cellstr(sprintf('%2.2f\t',P(:,1))), ' /',{' '}, cellstr(sprintf('%2.2f\t',P(:,2))) );
% fprintf('%s\t',str{:}); fprintf('%s\t\t','' ); fprintf('%s\t',str2{:}); fprintf('%s\n','' );
% add numbers on plot
% for j=1:length(P);   text(j, max(P(j,:))+max(E(j,:))+1, str{j}, 'horizontalalignment', 'center', 'fontsize', 8, 'parent', gca );   end



% Get a reduced Res table, just with BVCs/BCs in.
ResBVC = Res( any(Res.dataInd(:,1:2),2) & Res.isBVC & Res.area==area,   :   );

if 0
    % This is to get the highly reduced list for Edvard.
    fToKeep = {'cellID', 'rat', 'age', 'sessionID', 'SIAdSm', 'borderScoreAdSm'};
    f = Res.Properties.VariableNames;
    for ii=1:length(f)
        if ~strcmp( f{ii}, fToKeep )
            ResBVC.(f{ii}) = [];
        elseif size( Res.(f{ii}), 2 ) >= 2
            ResBVC.(f{ii}) = ResBVC.(f{ii})( :, 1:2 );
        end
    end
end

end

