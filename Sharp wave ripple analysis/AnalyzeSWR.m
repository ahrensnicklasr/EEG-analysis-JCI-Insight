function[]=AnalyzeSWR(fs, chNum, animalID)

%%Main funciton to loop through RHD files and run SWR_detect 
%% inputs include sampling rate,channel number, and animal and channel ID info

 d=dir('*.rhd');

disp(d);


z=1; %% index for all SWRs from and animal
for k=1:length(d)

 amp_data=command_read_Intan_RHD2000_file(d(k).name);
 cur_data=amp_data(chNum,:);

 %% transform data into single second rows
    a=1; 
    for i= 1:fs:(length(cur_data)-fs)
    sec_data(a,:)=cur_data(i:fs+i);
    a=a+1;
    end

%% calculate skew and rms of every second for the channel
    a=1; 
    for i=1:length(sec_data(:,1))
    skw(i)=skewness(sec_data(i,:));
    a=a+1;
    end

    a=1; 
    for i=1:length(sec_data(:,1))
    rtms(i)=rms(sec_data(i,:));
    a=a+1;
    end

%%find average and standard deviation of skw and rms
    stats(1)=mean(skw);
    stats(2)=std(skw);
    stats(3)=mean(rtms);
    stats(4)=std(rtms);
    
%% label bad channels (mean rms >200 or <30, or skw >0.4, provides 90% accuracy)   
    if stats(3)>200 || stats(3)<30 || stats(1)>0.4
        disp('bad channel');
        stats(5)=-1; %% flag for bad channels
    else 
        stats(5)=1; %% flag for good channel
    end
 
 if stats(5)>0   
 [clipped_data, zscore_p, SWR_index, SWR_voltages, zscore_SWR_p,  low_gamma_zscore_avg_SWR,  high_gamma_zscore_avg_SWR] = SWR_detect(cur_data, fs);
 summary(k,1)=length(clipped_data);
 summary(k,2)=length(SWR_index);
 all_clipped_data{z}(:)=clipped_data(:);
 all_zscore_p{z}(:,:)=zscore_p(:,:);
 
 
 all_SWR_voltages{z}(:,:)=SWR_voltages(:,:);
 final_SWR_voltages=zeros(1,2001);
 a=1; 
 for i=1:length(all_SWR_voltages) 
     if all_SWR_voltages{i}~=0
     final_SWR_voltages(a:(a+size(all_SWR_voltages{i},1)-1),:,:)=all_SWR_voltages{i}(:,:,:); 
     a=a+size(all_SWR_voltages{i},1);
     end
 end
 
 all_zscore_SWR_p{z}(:,:,:)=zscore_SWR_p(:,:,:);
 final_zscore_SWR_p=zeros(1,1024,11);
 a=1; 
 for i=1:length(all_zscore_SWR_p) 
     if all_zscore_SWR_p{i}~=0
     final_zscore_SWR_p(a:(a+size(all_zscore_SWR_p{i},1)-1),:,:)=all_zscore_SWR_p{i}(:,:,:); 
     a=a+size(all_zscore_SWR_p{i},1);
     end
 end
 
 all_low_gamma_zscore_avg_SWR{z}(:,:)=low_gamma_zscore_avg_SWR(:,:);
 final_low_gamma_zscore_avg_SWR=zeros(1,11);
 a=1; 
 for i=1:length(all_low_gamma_zscore_avg_SWR) 
     if all_low_gamma_zscore_avg_SWR{i}~=0
     final_low_gamma_zscore_avg_SWR(a:(a+size(all_low_gamma_zscore_avg_SWR{i},1)-1),:)=all_low_gamma_zscore_avg_SWR{i}(:,:); 
     a=a+size(all_low_gamma_zscore_avg_SWR{i},1);
     end
 end

 all_high_gamma_zscore_avg_SWR{z}(:,:)=high_gamma_zscore_avg_SWR(:,:);
 final_high_gamma_zscore_avg_SWR=zeros(1,11);
 a=1; 
 for i=1:length(all_high_gamma_zscore_avg_SWR) 
     if all_high_gamma_zscore_avg_SWR{i}~=0
     final_high_gamma_zscore_avg_SWR(a:(a+size(all_high_gamma_zscore_avg_SWR{i},1)-1),:)=all_low_gamma_zscore_avg_SWR{i}(:,:); 
     a=a+size(all_high_gamma_zscore_avg_SWR{i},1);
     end
 end

 z=z+1;
 end
 
 if stats(5)<0
 summary(k,1)=0/0;
 summary(k,2)=0/0;
 end
 
end
  
sum_filename=sprintf('SWRsummary_%s.mat',animalID);
save(sum_filename, 'summary');

clip_data_filname=sprintf('clipped_data_%s.mat',animalID);
save(clip_data_filname, 'all_clipped_data');
% 
% zscore_p_filname=sprintf('zscore_p_%s.mat',animalID);
% save(zscore_p_filname, 'all_zscore_p', '-v7.3');

% SWR_voltages_filename=sprintf('SWR_voltages_%s.mat', animalID);
% save(SWR_voltages_filename, 'all_SWR_voltages');

SWR_voltages_filename=sprintf('SWR_voltages_%s.mat', animalID);
save(SWR_voltages_filename, 'final_SWR_voltages');

% zscore_SWR_filename=sprintf('zscore_SWR_%s.mat', animalID);
% save(zscore_SWR_filename, 'all_zscore_SWR_p');

zscore_SWR_filename=sprintf('zscore_SWR_%s.mat', animalID);
save(zscore_SWR_filename, 'final_zscore_SWR_p');

% zscore_low_gamma_file =sprintf('low_gamma_zscore_SWR_%s.mat', animalID);
% save(zscore_low_gamma_file, 'all_low_gamma_zscore_avg_SWR');

zscore_low_gamma_file =sprintf('low_gamma_zscore_SWR_%s.mat', animalID);
save(zscore_low_gamma_file, 'final_low_gamma_zscore_avg_SWR');


% zscore_high_gamma_file =sprintf('high_gamma_zscore_SWR_%s.mat', animalID);
% save(zscore_high_gamma_file, 'all_high_gamma_zscore_avg_SWR');

zscore_high_gamma_file =sprintf('high_gamma_zscore_SWR_%s.mat', animalID);
save(zscore_high_gamma_file, 'final_high_gamma_zscore_avg_SWR');



end