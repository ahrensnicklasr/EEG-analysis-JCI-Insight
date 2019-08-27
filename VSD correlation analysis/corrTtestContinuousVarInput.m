function [pMat,tMat] = corrTtestContinuousVarInput(stackCtl,stackExp,xLenAvg,tLenAvg,debug)

if nargin==4,
    debug=0;
end

if debug,
    tic % start the timer
end

kernel = ones(xLenAvg,tLenAvg); % convolution kernel (for spatial and temporal averaging)

% preallocate
stackConvCtl=zeros(size(stackCtl));
stackConvExp=zeros(size(stackExp));


for i=1:size(stackCtl,3),
    stackConvCtl(:,:,i) = corrNanconvRepEdge(stackCtl(:,:,i),kernel);
end    

for i=1:size(stackExp,3),
    stackConvExp(:,:,i) = corrNanconvRepEdge(stackExp(:,:,i),kernel);
end

pMat=zeros(size(stackCtl(:,:,1))); % preallocate
tMat=zeros(size(stackCtl(:,:,1)));

for i=1:size(stackCtl,1),
    for j=1:size(stackCtl,2),
        valCtl= stackConvCtl(i,j,:); % get the convoluted values for the pixel (i,j) in the stack of rasters
        valExp= stackConvExp(i,j,:);
        
        valCtl(isnan(valCtl))=[]; % remove NaNs from the list
        valExp(isnan(valExp))=[];
        
        [~,pMat(i,j),~,stats] = ttest2(valCtl,valExp); % take stat across the stack of rasters
        tMat(i,j)=stats.tstat; % store the T statistic for each pixel
    end
end

if debug,

    toc % stop the timer
    
    figure(1212),set(gcf,'Name','rasterTtestContinuousVarInput','Position',[807 408 640 160],'PaperPositionMode','auto')
    imagesc(pMat)
    set(gca,'YDir','normal','CLim',[0 0.05]) % puts the first row of the raster at the bottom of the plot and sets colormap so everything with p>0.5 is one color
    box off
    
end

