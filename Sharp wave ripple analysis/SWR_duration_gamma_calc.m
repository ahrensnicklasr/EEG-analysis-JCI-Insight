function[SWR_duration] = SWR_duration_gamma_calc(SWR_voltages, low_gamma_zscore_SWR, fs)
%Function to take stack of SWR voltages and stack of low gamma during SWR
%and return matrix with duration of SWR and mean low gamma for 200ms
%following peak

f_SWR=[150, 250];
SWR_filt_voltages=zeros(size(SWR_voltages,1), length(SWR_voltages));
SWR_env=zeros(size(SWR_voltages,1), length(SWR_voltages));
SWR_env_zscore=zeros(size(SWR_voltages,1), length(SWR_voltages));
SWR_duration=zeros(size(SWR_voltages,1),2);

for i=1:size(SWR_voltages,1)
SWR_filt_voltages(i,:)=bandpass(SWR_voltages(i,:),f_SWR, fs);
[SWR_env(i,:), ~]=envelope(SWR_filt_voltages(i,:));
SWR_env_zscore(i,:)=zscore(SWR_env(i,:));
SWR_start_flag=-99;
SWR_end_flag=-99;
time_start=0;
time_end=NaN;
    for j=1000:1949
        if SWR_start_flag<0 && SWR_end_flag<0 && SWR_env_zscore(i,j)>3
            time_start=j;
            SWR_start_flag=1;
        end
        
        if SWR_start_flag>0 && SWR_end_flag<0 && SWR_env_zscore(i,j)<2
            
             if max(SWR_env_zscore(i,j:(j+50)))<2
                time_end=j;
                SWR_end_flag=1;
             end
          
        end
    end
    SWR_duration(i,1)=time_end-time_start;
    SWR_duration(i,2)=low_gamma_zscore_SWR(i,6);
end
SWR_duration(any(isnan(SWR_duration),2),:) = [];
end