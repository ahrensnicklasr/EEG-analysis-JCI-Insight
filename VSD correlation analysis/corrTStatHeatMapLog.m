function [Pmat,heatMap,colorbarForLegend,colorbarTickPositions,colorbarTickLabels] = corrTStatHeatMapLog(tMatReal,allTs,frameInterval_ms,tStim,t1,t2,axHandle,axHandleColormap)
% Given a 2D array of T-statistics TMATREAL and a list of 
% permutation-generated T-statistics ALLTS to compare against,
% TSTATHEATMAPLOG creates a pseudocolor heat map of P values.  
% P values <0.05 will be shaded purple.
%
% INPUTS:
% The inputs TMATREAL and ALLTS are required.  TMATREAL is the primary 
% output array from the function rasterPermTestCompare(). ALLTS is the 
% primary output array from the function rasterPermTestIterate().
%
%
% TSTATHEATMAPLOG(tMatReal,allTs,frameInterval_ms) compares the T
% statistic matrix TMATREAL (from the control vs. experimental test)
% to the emperical null distribution matrix ALLTS to compute and plot 
% P values for all sites in a raster plot.
%
% TSTATHEATMAPLOG(tMatReal,allTs,frameInterval_ms) plots the heatmap of 
% P values with the temporal axis incremented in time steps between movie 
% frames.  (For a recording obtained at 500 frames per second, 
% frameInterval_ms = 2 ms.)
%
% TSTATHEATMAPLOG(tMatReal,allTs,frameInterval_ms,tStim) adjusts the temporal axis so 
% that t = 0 ms corresponds to the time of stimulus delivery.
%
% TSTATHEATMAPLOG(tMatReal,allTs,frameInterval_ms,tStim,t1,t2) plots the 
% heatmap, showing only the samples during the time interval t1:t2.  Note 
% that t1 and t2 are in units of time, not sample numbers.
%
% TSTATHEATMAPLOG(tMatReal,allTs,frameInterval_ms,tStim,t1,t2,axHandle) 
% plots the heatmap as described above, in the axes defined by the handle 
% axHandle.
%
% It is also possible to omit some inputs to TSTATHEATMAPLOG, by providing
% empty arrays, [], for some inputs.  For example, to plot a heatmap in the
% axes defined by axis handle "axHandle" without providing tStim, t1, or
% t2, the following command could be used:
%
% tStatHeatMapLog(tMatReal,allTs,2,[],[],[],axHandle)
%
%
% OUTPUTS:
% A number of variables generated by tStatHeatMapLog can be returned, in
% the following order:
% [Pmat,heatMap,colorbarForLegend,colorbarTickPositions,colorbarTickLabels] = tStatHeatMapLog(...)
%
% These variables are defined as follows:
% 
% Pmat: P values for each spatiotemporal site (same size as TMATREAL).
%
% heatMap: the RGB color-indexed image of the heatMap.  this is the image 
% plotted by TSTATHEATMAPLOG.
%
% colorbarForLegend: RGB color-indexed colormap, of size [1000 3].  This is
% the colormap plotted by TSTATHEATMAPLOG.
%
% colorbarTickPositions: for the colorbar, used to set the positions of the
% X axis labels, e.g., set(gca,'XTick',colorbarTickPositions) .
%
% colorbarTickLabels: for the colorbar, used to set the values of the X axis
% labels, e.g., set(gca,'XTickLabel',colorbarTickPositions) .
%
%


ALPHA = 0.05;
P_MIN = 0.001;


% if any inputs aren't provided, set default values
if ~exist('frameInterval_ms','var') % if frameInterval_ms doesn't exist
    frameInterval_ms = 1; % define it as 1
elseif isempty(frameInterval_ms), % or if the variable exists but is empty
    frameInterval_ms = 1; % define it as 1
end

if ~exist('tStim','var') % if tStim doesn't exist
    tStim = 1; % define it as 1
elseif isempty(tStim), % or if the variable exists but is empty
    tStim = 1; % define it as 1
end

if ~exist('t1','var') % if t1 doesn't exist
    t1 = 1*frameInterval_ms; % define it as 1
elseif isempty(t1), % or if the variable exists but is empty
    t1 = 1*frameInterval_ms; % define it as 1
end

if ~exist('t2','var') % if t2 doesn't exist
    t2 = size(tMatReal,2)*frameInterval_ms; % define t2 according to the number of temporal samples in the data input
elseif isempty(t2), % or if the variable exists but is empty
    t2 = size(tMatReal,2)*frameInterval_ms; % define t2 according to the number of temporal samples in the data input
end

if ~exist('axHandle','var') % if axHandle doesn't exist
    figure % create a figure and axes
    set(gcf,'Position',[661 264 410 210])
    axHandle = axes;
elseif isempty(axHandle), % or if the variable exists but is empty
    figure % create a figure and axes
    set(gcf,'Position',[661 264 410 210])
    axHandle = axes;
end

if ~exist('axHandleColormap','var') % if axHandle doesn't exist
    figure % create a figure and axes
    set(gcf,'Position',[661 100 410 60])
    axHandleColormap = axes;
elseif isempty(axHandleColormap), % or if the variable exists but is empty
    figure % create a figure and axes
    set(gcf,'Position',[661 100 410 60])
    axHandleColormap = axes;
end
% end of input control 


t = ((t1):frameInterval_ms:t2-1) - tStim; % numbers for labeling the axis

allTs = sort(abs(allTs),3); % sort the array of Ts 

% This gives a heatmap that shows where the real T value falls on a % linear list of T values
% TmaxMat = squeeze(allTs(:,ind,end)); % use the maximum T value for each site.  assuming we have done a massive number of iterations, the real T statistic will be encompassed by the distribution of all Ts (and the real T statistic is expected to fall in the tail of the distribution if our ctl vs. exp grouping accounts for the majoritiy of the variability in the data)
% heatMap = abs(tMatReal(:,ind))./TmaxMat;

% Compute a p value for each site (corresponds to the area of the tail that
% falls beyond the real T statistic, for the distribution of T stats for
% each site
Pmat = zeros(size(tMatReal)); % preallocate
for i=1:size(Pmat,1),
    for j=1:size(Pmat,2), % for each spatiotemporal site,
        indTreal = find(allTs(i,j,:)>abs(tMatReal(i,j)),1,'first'); % get the index of the first T statistic that is greater than the real T statistic
        if isempty(indTreal),
            indTreal=length(allTs(i,j,:)); % if the real T statistic is bigger than any of the T values in the permuted data, set the found index to the maximum possible value
        end
        Pmat(i,j) = 1-indTreal/length(allTs(i,j,:)); % because the total area of the distribution is 1, the area of the tail (the p value) is the ratio of the number of points in the tail to the total number of points in the distribution
    end
end


heatMap = pMatToLogImgGrayPurple(Pmat,ALPHA,P_MIN); % call function (below) to convert matrix to RGB image

pMinExp = log10(P_MIN); % for colormap
colorbarLegendVect = 10.^(linspace(pMinExp,0,1000));
colorbarForLegend = pMatToLogImgGrayPurple(colorbarLegendVect,ALPHA,P_MIN); % make a colorbar


colorbarTickLabels = [0.001:.001:.01 .02:.01:.1 .2:.1:1]; % for labelling X axis of color legend

colorbarTickPositions = zeros(size(colorbarTickLabels)); % preallocate
for i=1:length(colorbarTickLabels),
    colorbarTickPositions(i) = find(colorbarLegendVect>=colorbarTickLabels(i),1)/1000; % find the tick position that corresponds to the current tick
end

colorbarTickLabels=num2cell(colorbarTickLabels);

indNoLabels = [2:9 11:13 15:18 20:27]; % indices where the 

for i=1:length(indNoLabels),
    colorbarTickLabels{indNoLabels(i)}=' ';
end



set(gca,'XTick',colorbarTickPositions)
set(gca,'XTickLabel',colorbarTickLabels)
box off
set(gca,'YTick',[])

heatMap = cat(3, heatMap(:,:,1), heatMap(:,:,2), heatMap(:,:,3) ); % add rows of zeros to mark anatomical boundaries in each channel

imagesc(t,[],heatMap,'Parent',axHandle,[0 1])

set(axHandle,'YDir','normal')
xlabel('time (ms)')



%% Y labels
szIm = size(Pmat);

boundaryRows=[0,0]; %because code for correlation maps without boundaries
tickInd = zeros(2,1); % preallocate.  there will be a tick mark for the first and last rows, plus a tick mark for the first row in each anatomical region
tickVal = zeros(2,1);

tickInd(1) = 1;
tickVal(1) = 1;

tickInd(end) = szIm(1);
tickVal(end) = szIm(1);

set(axHandle,'YTick',tickInd,'YTickLabel',tickVal) % display the y ticks and y tick labels on the figure

%% colorbar

xRange = [0 1];

imagesc(xRange,[],colorbarForLegend,'Parent',axHandleColormap)
set(axHandleColormap,'XTick',colorbarTickPositions)
set(axHandleColormap,'XTickLabel',colorbarTickLabels)
set(axHandleColormap,'YTick',[])
set(gcf,'Name','Heatmap color legend (logarithmic)')
box off


function [grayPurpleLogIm] = pMatToLogImgGrayPurple(Pmat,ALPHA,P_MIN)

logPmat = log10(1./Pmat); % take the inverse log of P.  this will expand the range of P values so we can see them better

% create color maps
grayCmap = linspace(165/255,255/255,256)';
grayCmap = [grayCmap grayCmap grayCmap];

purpleCmapR = linspace(55/255,155/255,256)';
purpleCmapG = linspace(30/255,150/255,256)';
purpleCmapB = linspace(110/255,235/255,256)';

purpleCmap = [purpleCmapR purpleCmapG purpleCmapB];

% for log plots, colormaps need to be flipped because small values become
% big logs, and vice versa
grayCmap=flipud(grayCmap);
purpleCmap=flipud(purpleCmap);


% to get better color depth in the purple (<ALPHA) range, we will
% create separate arrays for above-ALPHA and below-ALPHA 
PmatInsignif=logPmat; % copy the matrix into PmatInsignif
PmatInsignif(Pmat<ALPHA)=NaN; % NaN-out all the values that are significant
PmatSignif = logPmat;
PmatSignif(Pmat>=ALPHA)=NaN; % NaN-out all the insignificant values; PmatSignif now contains only significant values less than ALPHA; will be rendered purple

% convert each Pmat matrix to RGB values (result will be MxNx3)
PmatInsignifScale = PmatInsignif/log10(1/ALPHA); % stretch to use the full color scale
PimgGray = ind2rgb(gray2ind(PmatInsignifScale,256),grayCmap);

PmatSignifScale = (PmatSignif-log10(1/ALPHA))/(log10(1/P_MIN)-log10(1/ALPHA));
PmatSignifScale(isnan(PmatSignif))=0; % zero all sites that are supposed to be zero (subtraction on the line above caused zeros to become -log10(1/ALPHAs) )
PimgPurple=ind2rgb(gray2ind(PmatSignifScale,256),purpleCmap);

% combine the two images, one channel at a time

% first, set values to zero for sites that will not be colored with a given
% map (for example, sites with ALPHA<0.05 will be zeroed in the gray map)
PimgGrayR = PimgGray(:,:,1);
PimgGrayR(isnan(PmatInsignif)) = 0; % if a site was zero in PmatInsignif, make it zero in the image output
PimgGrayG = PimgGray(:,:,2);
PimgGrayG(isnan(PmatInsignif)) = 0; % same for G channel
PimgGrayB = PimgGray(:,:,3);
PimgGrayB(isnan(PmatInsignif)) = 0; % same for B channel

PimgPurpleR = PimgPurple(:,:,1);
PimgPurpleR(~isnan(PmatInsignif)) = 0; % if a site was zero in PmatSignif, make it zero in the image output
PimgPurpleG = PimgPurple(:,:,2);
PimgPurpleG(~isnan(PmatInsignif)) = 0; % same for G channel
PimgPurpleB = PimgPurple(:,:,3);
PimgPurpleB(~isnan(PmatInsignif)) = 0; % same for B channel

grayPurpleLogIm = cat(3,PimgGrayR+PimgPurpleR,PimgGrayG+PimgPurpleG,PimgGrayB+PimgPurpleB);












