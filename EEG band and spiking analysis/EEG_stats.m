function [bands spikes spikeflag summary]= EEG_stats(data,fs)
%bands is a matrix with n channels x m epochs x 5 bands (delta, theta,
%alpha, beta, gamma)
%spikes returns a matrix with n channels x m 5second epochs with the number
%of spikes per 5s epoch (these are true spikes screened for artifact)
%summary is a matrix which returns for each channel:1)mean delta 2)std
%delta 3)mean theta 4)std theta 5)mean alpha 6)std alpha 7)mean beta 8)std
%beta 9)mean gamma 10)std gamma 11)single channel spikes 12)multichannel
%spikes


%% read each channel of data
b=1; for j=1:length(data(:,1))
    cur_data=data(b,:);

%% transform data into single second rows
a=1; for i= 1:fs:(length(cur_data)-fs)
sec_data(a,:)=cur_data(i:fs+i);
a=a+1;
end

%% calculate skew and rms of every second for each channel
a=1; 
for i=1:length(sec_data(:,1))
skw(b,i)=skewness(sec_data(i,:));
a=a+1;
end

a=1; 
for i=1:length(sec_data(:,1))
rtms(b,i)=rms(sec_data(i,:));
a=a+1;
end

b=b+1;
end

%%find average and standard deviation of skw and rms
for i=1:length(data(:,1))
stats(i,1)=mean(skw(i,:));
stats(i,2)=std(skw(i,:));
stats(i,3)=mean(rtms(i,:));
stats(i,4)=std(rtms(i,:));
end

%% label bad channels (mean rms >200 or <30, or skw >0.4, provides 90% accuracy)
for i=1:length(data(:,1))
    if stats(i,3)>200 || stats(i,3)<30 || stats(i,1)>0.4
        disp('bad channel');
        disp(i);
        stats(i,5)=-1; %% flag for bad channels
    else 
        stats(i,5)=1; %% flag for good channel
    end
end

%% calculate z-score for rms
for i=1:length(rtms(:,1))
zscore_rtms(i,:)=zscore(rtms(i,:));
end

for i=1:length(data(:,1)) %% run half-max amplitude spike detector
    
    [spikes(i) spikeflag(i,:)]=mySpikes(data(i,:),5,fs);
    
    if stats(i,5)<=0 %% only keep good channels
        spikes(i)=0/0;
        spikeflag(i,:)=0/0;
    end
end

for i=1:length(data(:,1)) %% ignore spike if there is a false spike within 10 seconds of spike  
    for j=(fs*5):(length(spikeflag(i,:))-(fs*5))  %%ignore checking spikes in first and last 5 seconds because of end problems
        if spikeflag(i,j)==1 && min(spikeflag(i,(j-(fs*5)):(j+(fs*5))))<0
          spikes(i)=spikes(i)-1;
          spikeflag(i,j)=-2;
        end
    end
end
     
multi_lead_spikes=zeros(length(data(:,1)),1);
for i=1:length(data(:,1))
    for j=(1+(fs/10)):(length(spikeflag(i,:))-(fs/10))
        if spikeflag(i,j)==1 && sum(sum(spikeflag(:,(j-(fs/10)):(j+(fs/10)))))>=2 %% consider multilead spikes to be within 100ms of each other
        multi_lead_spikes(i)=multi_lead_spikes(i)+1;
        end
    end
end



%% calculate max rtms z-score for 5 second epochs to identify times of artifacts that will be excluded from band analysis
for i=1:length(zscore_rtms(:,1))
    a=1;
    for j=1:5:length(zscore_rtms(1,:))-5
        zscore_rtms_max(i,a)=max(zscore_rtms(i, j:j+5));
        a=a+1;
    end
end

%% take good channels and convert into 5 second epochs band analysis
b=1; for j=1:length(data(:,1))
    if stats(b,5)>0 %% only run on good channels
        cur_data=data(b,:);
        fs_epoch=fs*5;
        a=1; for i= 1:fs_epoch:(length(cur_data)-fs_epoch)
            epoch_data(a,:)=cur_data(i:fs_epoch+i);
            a=a+1;
        end
            
            
        for i=1:length(epoch_data(:,1))-1 
            
            if zscore_rtms_max(b,i)<3 %% only calculate spectra for 5s epochs with no movement artifacts (i.e. peak Z-score of RMS <3) 
                bands(b,i,1:5)=mySpectrum(epoch_data(i,:), fs);
            else
                bands(b,i,1:5)=0/0; %set flag for bad epoch
            end
           

        end
    else 
        bands(b,:,1:5)=0/0; %set flag for bad channel
    end
    b=b+1; 
end



for i=1:length(data(:,1))
    
    a=bands(i,:,1);
    b=bands(i,:,2);
    c=bands(i,:,3);
    d=bands(i,:,4);
    e=bands(i,:,5);
    f=spikes(i);
    g=multi_lead_spikes(i);
    
    summary(i,1)=mean(a(a>0)); %mean delta
    summary(i,2)=std(a(a>0)); %std delta
    summary(i,3)=mean(b(b>0)); %mean theta
    summary(i,4)=std(b(b>0)); %std theta
    summary(i,5)=mean(c(c>0)); %mean alpha
    summary(i,6)=std(c(c>0)); %std alpha
    summary(i,7)=mean(d(d>0)); %mean beta
    summary(i,8)=std(d(d>0)); %std beta
    summary(i,9)=mean(e(e>0)); %mean beta
    summary(i,10)=std(e(e>0)); %std beta
    summary(i,11)=sum(f(f>0)); %total number of spikes
    summary(i,12)=sum(g(g>0)); %total number of shared spikes
    if sum(isnan(summary(i,:)))>0
        summary(i,11)=NaN;
        summary(i,12)=NaN;
    end
    end



end