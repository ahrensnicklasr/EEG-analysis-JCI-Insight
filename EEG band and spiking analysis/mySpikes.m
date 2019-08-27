
function[spikes spikeflag]=mySpikes(data,Zspike, fs)
%Spike detector using width at half maximum approach
%Called by EEG_stats
%Inputs are: 
    %data= single channel EEG data, 
    %Zspike= Zscore threshold to call something a spike
    %fs= sampling frequency

%spikes are a minimum of 200ms apart
    %Outputs are:
    %spikeflag 1=spike, 0=no spike, -1=artifact
    %spikes= total number of spikes
    
[b,a]=butter(6,((70/fs)*2)); %Lowpass filter at 70Hz
data=filter(b,a,data);
spikeflag=zeros(1,length(data));
spikes=0;
zData=0; 
zData=zscore(data); 
spikeThreshold=Zspike*std(data); %threshold to call a big deflection an event
[zDataPks zDataLocs]=findpeaks(zData,'MinPeakHeight', Zspike, 'MinPeakDistance', 0.2*fs); %find high amplitude deflections, separated by at lest 200ms

for i=1:length(zDataLocs)  
    a=0;
    j=zDataPks(i);
    
    halfWidthEnd=1;
    halfWidthBegin=1;
    
    while(j>(zDataPks(i)/2) && halfWidthEnd<length(zDataLocs))
        j=data(zDataLocs(i)+a);
        halfWidthEnd=zDataLocs(i)+a;
        a=a+1;
    end
    
    b=0;
    j=zDataPks(i);
    while(j>(zDataPks(i)/2) && halfWidthBegin>0)
        j=data(zDataLocs(i)-b);
        halfWidthBegin=zDataLocs(i)-b;
        b=b+1;
    end
    
    halfWidth(i,1)=(halfWidthEnd-halfWidthBegin); %Calculate maximum half width
  
    if halfWidth(i,1)>(0.005*fs) && halfWidth(i,1)<=(0.2*fs) %Spike max half width should be between 5-200 ms
        spikes=spikes+1;
        spikeflag(zDataLocs(i))=1;
    end
     if halfWidth(i,1)<=(0.005*fs) || halfWidth(i,1)>(0.2*fs) %If peak is not true spike set flag for bad spike
        spikeflag(zDataLocs(i))=-1;
    end
    
   
   
end

end