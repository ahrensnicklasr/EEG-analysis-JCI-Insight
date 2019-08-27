function [] = corrMatrixPermTestIterate(TempStackCtl,TempStackExp,numIterations,conditionID,debug)
%
% corrPERMTESTITERATE computes t-statitics for continuous small regions 
% that span the entire data field, for random reshufflings of the input 
% data sets STACKCTL and STACKEXP.
% 
% corrPERMTESTITERATE will use the distribution of T statistics for 
% each column to determine the proper cutoff T statistic that indicates
% significance in that column.
% 
% INPUTS:
%
% STACKCTL and STACKEXP are 3D stacks of raster data that contain control
% and experimental data, respectively.
%
% NUMITERATIONS is the number of reshufflings to perform (500 is a typical
% value).
%
% DEBUG, when set to 1, will cause RASTERBOOTSTRAPITERATE to plot the image
% of T-statistics created during each iteration.
%
%
% OUTPUTS:
%
% ALLTS is a 2D matrix containing all the T statistics computed.  ALLTS is 
% of the form [iterationNumber AllTstatisticsForThatIteration].
%
% TCUTOFF is the T value that corresponds to alpha=0.05 for the input data.
% This value determines which regions are and are not significantly
% different.  TCUTOFF will be used as input to RASTERBOOTSTRAPCOMPARE.
%

tic % start timer

if nargin==3,
    debug=0;
end

%% Parameters for permutation testing, used by rasterPermTestIterate() and rasterPermTestCompare().
N_X_SAM = 3; % number of spatial samples to average for measurement at each raster pixel.  this value must be odd.
N_T_SAM = 3; % number of temporal samples to average for measurement at each raster pixel.  this value must be odd.



N_X_SAM=N_X_SAM; N_T_SAM=N_T_SAM; % necessary for compatibility with parfor loop


stackCtl=permute(TempStackCtl, [3,2,1]);
stackExp=permute(TempStackExp, [3,2,1]);

% get the number of rasters in each stack
nCtl = size(stackCtl,3);
nExp = size(stackExp,3);

allTs = zeros(size(stackCtl,1),size(stackCtl,2),numIterations); % preallocate to store a *2D array* of T values from each iteration

parfor i=1:numIterations,
    combinedStack = cat(3,stackCtl,stackExp); % combine the Ctl and Exp stacks into one stack
    
    indShuffled = randperm(nCtl+nExp); % shuffled indices of all rasters
    randIndCtl=indShuffled(1:nCtl); % take only the first nCtl rasters in each stack, so a group of rasters is chosen of the same cardinality as stackCtl
    randIndExp=indShuffled(nCtl+1:end); % this will give a group of rasters of the same cardinality as stackExp
    
    permGroup1 = combinedStack(:,:,randIndCtl); % use half of the rasters in one group (the first half of the stack)
    permGroup2 = combinedStack(:,:,randIndExp); % put the other half of the rasters in another group (the second half of the stack)
    
    [~,tMat] = corrTtestContinuousVarInput(permGroup1,permGroup2,N_X_SAM,N_T_SAM); % compute the t-statistic to compare the two groups at each site
    
    allTs(:,:,i)=tMat; % save all the T values from the current iteration as a 2D matrix
    
    if ~rem(i,10),
        disp(i)
    end
    
end

tElapsed = toc; % save elapsed time
disp(['corrPermTestIterate: Elapsed time is ' num2str(tElapsed/60) ' minutes.']) % display elapsed time

[tMatReal, ~] = corrPermTestCompare(stackCtl,stackExp,allTs,0.05);
corrTStatHeatMapLog(tMatReal, allTs, [], []);

tMatReal_filename=sprintf('tMatReal_%s.mat', conditionID);
save(tMatReal_filename, 'tMatReal');

allTs_filename=sprintf('allTs_%s.mat', conditionID);
save(allTs_filename, 'allTs');

if debug,
    % display each 2D image in ALLTS, waiting for the user to press a button between each image
    figure(987)
    set(gcf,'Name','tMat from current iteration')
    for i=1:numIterations, % for each slice
        figure(987)
        imagesc(abs(allTs(:,:,i))) % use ABS so that the image is a "heat map" of the size of T values (would be harder to interpret with negative T values)
        axis xy
        display(['max T for iteration ' num2str(i) ' = ' num2str(max(max(abs(allTs(:,:,i)))))]) % print the highest T value in the current image
        pause
    end
end




