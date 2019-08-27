function[clipped_data, zscore_p, SWR_index, SWR_voltages, zscore_SWR_p,  low_gamma_zscore_avg_SWR,  high_gamma_zscore_avg_SWR] = SWR_detect(data, fs)
%Function to take data, filter, clip out movement artifact, indentify theta and non-theta periods, then detect
%SWR during non-theta periods and save. Adapted from publication Iccarino et
%al. PMID:27929004.

%% filter raw data and take RMS
fpass=[1 300]; %frequency for bandpass filter of all data
filt_data= bandpass(data,fpass,fs); %filter raw data
rms_filt_data= rms(filt_data);

%% divide filtered data into 5sec epochs and take RMS for each epoch
a=1; 
for i= 1:(fs*5):(length(filt_data)-(fs*5))
epoch_data(a,:)=filt_data(i:(fs*5)+i);
rms_epoch_data(a)=rms(epoch_data(a,:));
a=a+1;
end

%% calculate z-score of rms voltage for each 5sec epoch
zscore_rms_epoch=zscore(rms_epoch_data(:);

%% set epochs with rms zscore > 2 as artifact
for b=1:length(zscore_rms_epoch) 
    if zscore_rms_epoch(b)>2 
        epoch_data(b,:)=-9999;
    end
end
%% define new clipped, filtered, data array with movement artifact 10s epochs removed, use this dataset for remainder of function
clipped_data= reshape(epoch_data',[], 1);
clipped_data=clipped_data(clipped_data>-9998);

%%calculate spectrogram of entire clipped data recording to use for
%%normalization of SWR periods
[p,f,t]=pspectrum(clipped_data,2000,'spectrogram','FrequencyLimits',[1 300], 'TimeResolution', 0.1);

for i=1:length(f)
zscore_p(i,:)=zscore(p(i,:));
end

%% filter for each freq band for theta and SWR analysis 
f_delta= [1 4];
f_theta= [4 12];
f_beta=[12 30];
f_SWR=[150 250]; %frequency range defined for sharp wave ripples

delta_filt_data=bandpass(clipped_data,f_delta,fs);
theta_filt_data=bandpass(clipped_data,f_theta,fs);
beta_filt_data=bandpass(clipped_data,f_beta,fs);
SWR_filt_data=bandpass(clipped_data,f_SWR,fs);

%calculate envelope of data filtered to each frequency band
[delta_upper_env, delta_lower_env]=envelope(delta_filt_data);
[theta_upper_env, theta_lower_env]=envelope(theta_filt_data);
[beta_upper_env, beta_lower_env]=envelope(beta_filt_data);
[SWR_upper_env, SWR_lower_env]=envelope(SWR_filt_data);

%calculate theta ratios using filtered data
for i=1:length(clipped_data)
theta_ratio(i)=(theta_upper_env(i)/(beta_upper_env(i)+delta_upper_env(i)));
end

%calculate z-scores for the theta_ratio and SWR 
z_theta_ratio= zscore(theta_ratio);
z_SWR_upper_env= zscore(SWR_upper_env);

%% define SWR timepoints if theta z-score <-0.5, and SWR z-score >4
SWR_flag=zeros(1, length(clipped_data));
for i=1:length(clipped_data)
    if z_theta_ratio(i) <=-0.5 && z_SWR_upper_env(i) >4
        SWR_flag(i)=1;  
    else
        SWR_flag(i)=0;
    end
end

%% calculate a sum of spike flag over 30ms moving window
for i=((fs/1000)*15)+1:(length(SWR_flag)-((fs/1000)*15))       
    y(i)= sum(SWR_flag(i-((fs/1000)*15):i+((fs/1000)*15)));
end

%% Find time points of max SWR z-score where 10ms out of 30ms window with positive SWR flag (i.e. if theta z-score <-0.5, and SWR z-score >4), max true SWR z-score peaks must be separated by 5 sec
a=1;
i=fs*1;
while i<(length(SWR_flag)-(fs*11))    
    if y(i)>(fs/1000*10)
        [psor, lsor]=findpeaks(z_SWR_upper_env(i:i+(10*fs)));
        SWR_index(a)= lsor(1)+i;
        i=i+(5*fs);
        a=a+1;       
    else
        i=i+1;
    end
    
end


%% if a SWR is detected, find time, and save voltages, SWR filt voltages, power z-scores in 1sec surrounding peak
check_SWR_index= exist('SWR_index');
if check_SWR_index==1
    time_SWR_index(:)=SWR_index./fs;
    [q, f_index_gamma_start]=min(abs(20-f)); 
    [q, f_index_gamma_stop]= min(abs(50-f));
    [q, f_index_high_gamma_stop]= min(abs(80-f));
    for i=1:length(SWR_index)
     SWR_voltages(i,:)=clipped_data(SWR_index(i)-fs/2:SWR_index(i)+fs/2);
     SWR_filt_voltages(i,:)=SWR_filt_data(SWR_index(i)-fs/2:SWR_index(i)+fs/2);
     
     %%find time points surrounding SWR_index time in z-score scaled power spectrum
     [m,index] = min(abs(t-time_SWR_index(i)));
     zscore_SWR_p(i,:,:)=zscore_p(:,index-5:index+5); %save 10, 100ms duration z-score bins around peak
     low_gamma_SWR_zscore_p(i,:,:)= zscore_SWR_p(i,f_index_gamma_start:f_index_gamma_stop,:);
     low_gamma_zscore_avg_SWR(i,:)=mean(low_gamma_SWR_zscore_p(i,:,:),2);
     high_gamma_SWR_zscore_p(i,:,:)= zscore_SWR_p(i,f_index_gamma_stop:f_index_high_gamma_stop,:);
     high_gamma_zscore_avg_SWR(i,:)=mean(high_gamma_SWR_zscore_p(i,:,:),2);
    end
   
else 
    SWR_index=0;
    time_SWR_index=0;
    SWR_voltages=0;
    SWR_filt_voltages=0;
    zscore_SWR_p=0;
    low_gamma_SWR_zscore_p= 0;
    low_gamma_zscore_avg_SWR=0;
    high_gamma_SWR_zscore_p= 0;
    high_gamma_zscore_avg_SWR=0;
end



end