function [ResSim] = generateAllBVCResponses_Sim(Res,ModMaps)
% Generate simulated BVCs, using real path data and random poisson spikes,
% where spike probability is determined by a particular BVC tuning curve.
%
% As input, requires the Res table and model maps (for fitting), but also 
% requires raw data in SCAn is loaded in workspace.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter list
prms.nSimBVCPerTrial          = 500;   % How many BVCs are simulated per trial (i.e. per set of position data).
prms.trialsToUse              = 1;   % Which trials to use, for each cell
prms.ageBins                  = [16 18; 19 21; 22 25; 40 40];  % Need age bins to match mean rates by age group.
% Scaling data to fit standard square %
prms.minOccForEdge            = 50;
prms.boxExtent                = {250, 280};  %% There are two values, one for the 62.5cm regular box, and the other for the 70 cm Bjerknes box. When set to 250 and 280, this imples a ppm of 400.
% BVC tuning curve parameters
% prms.d_list                   = 1:13 * 10;                % Only running in squares for now.
prms.phi_list                 = deg2rad( 0:6:354 );  % Angular tunings. In radians.
% prms.sigZeros                 = [6.2   12.2   20.2   30.2];    % (2.4928 : 1 : 6.4928) .^ 2; % sqrt(sig0) varies linearly and sqrt(12.2)=3.4928
prms.kappa                    = 25; 
prms.beta                     = 183;
prms.binSize                  = 0.25;     % Bin size in cm. NOTE! this is set to 1 cam pixel, as simulation works at this resolution.
% Definitions of environment size/shape %
prms.mapSize                  = [252 252];  % All data is squares scaled to 250x250 pix, then add 10 pix boundary for wall (analogous to 1 bin for wall in binned real data).
prms.barrStr                  = [1      1     1    252    0; ...  % Barriers are defined as being in centre of 
                                 1     252   252   252    0; ...  % wall bins (which are themselves 10x10 pix)
                                 252   252   252    1     0; ...  % again analogous to real binned data.
                                 252    1     1     1     0];
% Load in mat file with existing BVC props, for those characteristics which 
% should match real data population.
S = load('BVCPropsForSim.mat');   BVCPropsForSim = S.BVCPropsForSim;
% NB! Saved d props are in cm. For simulation, need to convert to cam pix
BVCPropsForSim.dTuning.data = (BVCPropsForSim.dTuning.data .* 4); 

% Convert the d and phi actual bins into an index (for rebinning corr outputs) %
[~,~,dBinInd]   = unique( ModMaps.dKey );
[~,~,phiBinInd] = unique( ModMaps.phiKey );
[~,~,szBinInd]  = unique( ModMaps.szKey );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce the rows in Res so that we've got just one row per SUB dataset 
% (i.e. remove 'duplicate' cells recorded from the same trial).
Res            = Res( Res.area==1, : );
Res            = Res( Res.age>15,  : );
[dsList,dsInd] = unique( Res.dataset, 'stable' );
Res            = Res( dsInd, : );
% For PARFOR, can't use table directly so need to extract some variables, also
% create containers for others, etc.
ageByDS        = Res.age;
XY             = cell(height(Res),prms.trialsToUse);
bvcMapSim      = cell(height(Res),prms.trialsToUse,prms.nSimBVCPerTrial);
simBVCProps    = nan(height(Res), 4, prms.nSimBVCPerTrial);  % Props are the same across trials (if multiple trials run)
fitBVCProps    = nan(height(Res), 4, prms.nSimBVCPerTrial);  % so 2nd dim used to store multiple values.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Gather the necessary data. Get a list of trials from Res, and add prms.phi_list
%     raw X and Y pixel-resolution data from scan.
for itDS = 1:length(dsList)
    % For this dataset, get from base workspace, then apply standard scaling.
    data      = evalin('base',dsList{itDS});
    data      = devSub_findEdgesAndScalePath( data, prms.minOccForEdge, prms.boxExtent{ 1 } );  % boxExtent always 250, as only running this for SUB data
    % Loop over cells/trials, and assign raw XY to cell array, such that indexing matches Res table.
    for itTr = 1:length( prms.trialsToUse ) 
        trIndInScan       = Res.trialInd(itDS,prms.trialsToUse(itTr));
        if isnan(trIndInScan);   continue;    end   % If a particular trial type e.g. CCE has not been run in this expt.
        XY{ itDS, itTr }  =  [data.trials( trIndInScan ).x, data.trials( trIndInScan ).y];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Loop over trial data and simulate BVCs
itNSVect = 1:prms.nSimBVCPerTrial;  % Need to construct iteration vectors for use inside the parfor loop, outside of it.
parfor itDS=1:size(XY,1)
    
    % Get distributions of tuning properties from which to sample: for phi this is flat,
    % for d, SZ, and rate, these are drawn from actual data distributions.
    ageBinForDS    = find( ageByDS(itDS)>=prms.ageBins(:,1) & ageByDS(itDS)<=prms.ageBins(:,2) );
    distDistForAge = BVCPropsForSim.dTuning.data( BVCPropsForSim.dTuning.ageKey==ageBinForDS );
    phiDistForAge  = BVCPropsForSim.phiTuning.data( BVCPropsForSim.phiTuning.ageKey==ageBinForDS );
    SZDistForAge   = BVCPropsForSim.sigZero.data( BVCPropsForSim.sigZero.ageKey==ageBinForDS );
    MRDistForAge   = BVCPropsForSim.meanRate.data( BVCPropsForSim.meanRate.ageKey==ageBinForDS );
    nTuneParams    = [length(distDistForAge), length(phiDistForAge), length(SZDistForAge), length(MRDistForAge)];
    
    for itNS=itNSVect

        % Randomly sample tuning characteristics for this BVC:
        dist     = distDistForAge( randi(nTuneParams(1)) ); 
        phi      = phiDistForAge( randi(nTuneParams(2)) );
        sz       = SZDistForAge( randi(nTuneParams(3)) );
        meanRate = MRDistForAge( randi(nTuneParams(4)) );
        tuning   = [prms.beta, sz, prms.kappa];    
        
        % Run simulation and fit (potentially on more than one trial, to check stability etc.
        for itTr=1 % length(prms.trialsToUse)
                  
            % Run simulation
            bvcMapSim{itDS,itTr,itNS} = makeBVCMap_v2( dist, phi, prms.mapSize, prms.barrStr, [], prms.binSize, tuning, XY{itDS}, meanRate );
            simBVCProps(itDS,:,itNS)  = [dist, phi, sz, meanRate];

            % Run the simulated BVC through the fitting.
            BVCModMaps   = ModMaps.modelMaps{ strcmp( {ModMaps.mapDefs(1:end).env}, 'hp' ) };
            rateMap      = padarray( bvcMapSim{itDS,itTr,itNS}, [1 1], NaN );
            % Get the pearson's r by hand so as to vectorise, this is the equation:
            A     = rateMap - nanmean(rateMap(:));                                      % Subtract samples sample mean. Dims of A here are 1,2
            B     = bsxfun(@minus, BVCModMaps, nanmean( nanmean(BVCModMaps,1), 2 ) );   % Dims of B are 1,2,3, mean for subtraction is calculated along third dim.
            AB    = bsxfun(@times, A, B);
            sumAB =  nansum(nansum(AB,1),2);                    % This is a 1x1xnMod vector, one sum(AB) for each model.
            sqrAA = sqrt( nansum(nansum(A.^2,1),2) );           % This is a 1x1 integer
            sqrBB = sqrt( nansum(nansum(B.^2,1),2) );           % This is a 1x1xnMod vector, one sqrt(B.^2) for each model.
            R     = sumAB ./ bsxfun(@times, sqrAA, sqrBB);      % This is a 1x1xnMod vector, one pearsons-r for each model.
            RArr  = accumarray( [dBinInd, phiBinInd, szBinInd], squeeze(R), [max(dBinInd), max(phiBinInd), max(szBinInd)], @nansum, nan );   % Reshape the r-vals into a ( 1:nPhiBin, 1:ndBin, 1:szBin ) array. (note that the 'nansum' isn't actually summing anything, there's just one value per bin, it's just a way of assigning values robust to NaNs.
            % Get best fit.
            tempProps                                                = nan(1,4);
            [tempProps(4), tempProps(1), tempProps(2), tempProps(3)] = nanmax2D(RArr);
            fitBVCProps(itDS,:,itNS)                                 = tempProps;

%             figure;  imagesc(bvcMapSim{itDS,itTr});
%             title(   sprintf('d=%d, phi=%d, SZ=%3.2f',dist,round(rad2deg(phi)),sz)  );


        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile output into a table
ResSim = table( dsList, ageByDS, fitBVCProps, simBVCProps, bvcMapSim, 'VariableNames', {'dataset','age','fitProps','simProps','rateMap'});













