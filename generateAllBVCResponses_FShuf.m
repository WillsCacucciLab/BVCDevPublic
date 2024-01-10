function [ResBVC] = generateAllBVCResponses_FShuf(Res,ModMaps,varargin)
% Batch function for generating BVC responses for all data, for randomised data, using .
%
%       [Res] = generateAllBVCResponses_FShuf(cellIDList,ModMaps,DataIn,varargin)
%
% Inputs: Res      - Existing results from getBasicSpatScores.
%         ModMaps  - canonical BVC maps structure
%
% Outputs: Res     - Same results table as input, with added fields:
%                     - BVCResp_fShuf, correlation of neural rate map with all model BVC fields. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Parameters %
% Stuff on how this function itself works.
prms.showTimer       = 1;
prms.makeCircRespMap = 0;
prms.boxExtent       = Res.Properties.UserData.boxExtent;
prms.BVCrBin         = Res.Properties.UserData.BVCrBin;
prms.nShuf           = 200;

% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for itCl=1:2:length(varargin);   prms.(varargin{itCl}) = varargin{itCl+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for itCl=1:length(f);   prms.(f{itCl}) = s.(f{itCl});   end                  %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup %%%


% Need to make copies of vars that cannot reference directly inside PARFOR loop.
cellIDList              = Res.cellID; 
envTypeList             = Res.envType;
ratList                 = Res.rat;
isMoserData             = Res.area==2;
rateMaps_unSm           = Res.rateMap_UnSm;   % Used for field shuffle
rateMaps_pos            = Res.posMap;    %   "
maxNTrial          = 5;
trialsToUse        = [1 2];


% Initialise cell arrays to contain results
cellID               = cell( size(cellIDList) );
ShufBVCRespMax_VFS   = cell( size(Res.rateMap)); 


WaitMessage = parfor_wait( length(cellIDList) );
tic
parfor itCl=1:length(cellIDList)
    WaitMessage.Send;
    if rem(itCl,100)==0
        disp(itCl);
    end

    BVCMaxForCell = cell( 1, maxNTrial );

    isEnvSqFlag     = ~strcmp( envTypeList(itCl, :), 'CCE');

    mapsForCl_unSm  = rateMaps_unSm(itCl,:);
    mapsForCl_pos   = rateMaps_pos(itCl,:);

    for itTr=1:length(trialsToUse)
        
        if isnan(trialsToUse(itTr));   continue;   end
        BVCMax_fShuf = nan(prms.nShuf,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Before either analysis, need get the environment size and type. Need to take account of:
        % - Moser square boxes larger than london boxes
        % - some London rats odd barriers
        % - is the data transformed square->circle and vice versa?
        % i) Size of env (max map dimension).
        if ~isEnvSqFlag(  trialsToUse(itTr)  )
            mapSize = prms.boxExtent(3)/prms.BVCrBin;   % All circles (either 'real' or transformed squares) are one standard size by design.
        elseif isMoserData(itCl)
            mapSize = prms.boxExtent(2)/prms.BVCrBin;   % All Moser squares are 280 pix side length
        else
            mapSize = prms.boxExtent(1)/prms.BVCrBin;   % All London squares are 250 pix side length
        end
        if isempty(mapsForCl_unSm{ trialsToUse(itTr) });  continue;  end

        Map_unSm   = mapsForCl_unSm{ trialsToUse(itTr) }(1:mapSize, 1:mapSize);
        Map_pos    = mapsForCl_pos{ trialsToUse(itTr) }(1:mapSize, 1:mapSize);

        if sum(Map_unSm,"all","omitnan") == 0;  continue;  end  % These shouldn't be in the dataset? But they are for some reason, and crash the voronoi func


        % ii) Get the right set of BVC models depending on experiment + use of shape transform
        if ~isEnvSqFlag(  trialsToUse(itTr)  )
            % i. First, if shape flag says circle, use circle
            EnvType = 'CCE'; 
        elseif strcmp( envTypeList(itCl,trialsToUse(itTr)), 'CCE' )
            % ii. If env label in raw data says 'CCE', but shape flag is square, data has been transformed - use baseline square.
            EnvType = 'hp';
        elseif strncmp( envTypeList(itCl,trialsToUse(itTr)), 'barrier', 7 )   &&   any( strcmp({'123', '124', '137', '138'}, ratList{itCl}) )
            % iii. Some London rats were run with shorter barriers, but these are not labeled specifcally in raw data.
            EnvType = [ envTypeList{itCl,trialsToUse(itTr)}, ' short' ];
        else
            % iv. In all other cases, we can trust raw data label to identify correct model map.
            EnvType = envTypeList{itCl,trialsToUse(itTr)};
        end
        BVCModMaps   = ModMaps.modelMaps{ strcmp( {ModMaps.mapDefs(1:end).env}, EnvType ) };
        

        % Make the shuffled maps% 
        mapsShuf     = voronoiFieldShuffle( Map_unSm, Map_pos, prms.nShuf, cellIDList{itCl}, itTr );
        
        for itSh = 1:prms.nShuf
            
            rateMap =  padarray( mapsShuf(:,:,itSh), [1 1], NaN );

            % Get the pearson's r by hand so as to vectorise, this is the equation:
            % a = a - mean2(a);     b = b - mean2(b);
            % r = sum(sum(  a.*b  ))        /       (  sqrt( sum(sum(  a.*a  ))  *   sqrt( sum(sum(   b.*b   )) )    );
            A  = rateMap - nanmean(rateMap(:));                                      % Subtract samples sample mean. Dims of A here are 1,2
            B  = bsxfun(@minus, BVCModMaps, nanmean( nanmean(BVCModMaps,1), 2 ) );   % Dims of B are 1,2,3, mean for subtraction is calculated along third dim.
            AB = bsxfun(@times, A, B);
            sumAB =  nansum(nansum(AB,1),2);                    % This is a 1x1xnMod vector, one sum(AB) for each model.
            sqrAA = sqrt( nansum(nansum(A.^2,1),2) );           % This is a 1x1 integer
            sqrBB = sqrt( nansum(nansum(B.^2,1),2) );           % This is a 1x1xnMod vector, one sqrt(B.^2) for each model.
            R     = sumAB ./ bsxfun(@times, sqrAA, sqrBB);      % This is a 1x1xnMod vector, one pearsons-r for each model.
            BVCMax_fShuf(itSh) = max(R,[],"all","omitnan");

        end

        BVCMaxForCell{ 1, trialsToUse(itTr) } = BVCMax_fShuf; 

        % END for by-trial loop.
    end        


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) Assign results for this dataset to output table % 

    cellID( itCl, 1 )              = cellIDList(itCl);
    ShufBVCRespMax_VFS( itCl, : )  = BVCMaxForCell; 
        
end
WaitMessage.Destroy;
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Merge BVC response results with existing Res table.
ResBVC = table( cellID, ShufBVCRespMax_VFS );
save( 'vs_Res_temp.mat', 'ResBVC' );

Res    = join( Res, ResBVC, "Keys", 'cellID', 'rightvariables', 'ShufBVCRespMax_VFS');










