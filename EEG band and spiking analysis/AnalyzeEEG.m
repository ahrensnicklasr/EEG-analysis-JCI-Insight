function[]=AnalyzeEEG(fs, numCages, cg1st, cg1end, cg2st, cg2end, cg3st, cg3end, cg4st, cg4end)

%%Main funciton to loop through RHD files and run EEG_stats which calls
%%spike detector (mySpikes) and spectrum analysis (mySpectrum)
%% inputs include sampling rate, number of cages, channel # for start of cage 1, channel# for end of cage 1, channel # for start of cage2, channel# for end of cage 2, channel# for start of cage 3, channel# for end of cage 3, channel# for start of cage4, channel# for end of cage4

 d=dir('*.rhd');

disp(d)

a=1;
for i=1:length(d);
 amp_data=command_read_Intan_RHD2000_file(d(i).name);
  
    if numCages>0
        [cage1bands{i} cage1spikes{i} cage1spikeflag{i} cage1summary{i}]=EEG_stats(amp_data(cg1st:cg1end,:), fs);
    end

    if numCages>1
       [cage2bands{i} cage2spikes{i} cage2spikeflag{i} cage2summary{i}]=EEG_stats(amp_data(cg2st:cg2end,:), fs);
    end
    

    if numCages>2
        [cage3bands{i} cage3spikes{i} cage3spikeflag{i} cage3summary{i}]=EEG_stats(amp_data(cg3st:cg3end,:), fs);
    end
    
    if numCages>3
        [cage4bands{i} cage4spikes{i} cage4spikeflag{i} cage4summary{i}]=EEG_stats(amp_data(cg4st:cg4end,:), fs);
    end
    
    
    
    if(numCages>0)
        save('cage1bands.mat', 'cage1bands');
        save('cage1spikeflag.mat', 'cage1spikeflag', '-v7.3');
        save('cage1summary.mat', 'cage1summary');
    end
    if(numCages>1)
        save('cage2bands.mat', 'cage2bands');
        save('cage2spikeflag.mat', 'cage2spikeflag', '-v7.3');
        save('cage2summary.mat', 'cage2summary');
    end
    if(numCages>2)
        save('cage3bands.mat', 'cage3bands');
        save('cage3spikeflag.mat', 'cage3spikeflag', '-v7.3');
        save('cage3summary.mat', 'cage3summary');
    end
    if(numCages>3)
        save('cage4bands.mat', 'cage4bands');
        save('cage4spikeflag.mat', 'cage4spikeflag', '-v7.3');
        save('cage4summary.mat', 'cage4summary');
    end
    a=a+1;
end

end