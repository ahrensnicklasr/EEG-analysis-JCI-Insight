function[]=AnalyzeVSDIcorr(stretchStack, conditionID)

%%Main funciton to 1) take a stretched VSDI raster stack, and condition label, 2) calculate
%%correlation coefficient stack in time and space, 3) collapse to
%%average correlation matrix 4) calculate and save regional avg r values
%%(Hilus to Hilus, Hilus to CA3, Hilus to CA1, CA3 to CA3, CA3 to CA1, CA1
%%to CA1)

for i=1:length(stretchStack)
    corr_time(i,:,:)=corr(stretchStack{i});
    corr_space(i,:,:)=corr(stretchStack{i}');
    region_avg(i,1)=mean(mean(corr_space(i,1:5,1:5))); %hilus-hilus corr
    region_avg(i,2)=mean(mean(corr_space(i,1:5,6:18))); %hilus-CA3 corr
    region_avg(i,3)=mean(mean(corr_space(i,1:5,19:end))); %hilus-CA1 corr
    region_avg(i,4)=mean(mean(corr_space(i,6:18, 6:18))); %CA3-CA3 corr
    region_avg(i,5)=mean(mean(corr_space(i,6:18, 19:end))); %CA3-CA1 corr
    region_avg(i,6)=mean(mean(corr_space(i,19:end, 19:end))); %CA1-CA1 corr
end
mean_corr_space(:,:)=mean(corr_space, 1);
mean_corr_time(:,:)=mean(corr_time, 1);

time_corr_stack_filename=sprintf('time_corr_stack_%s.mat', conditionID);
save(time_corr_stack_filename, 'corr_time');

space_corr_stack_filename=sprintf('space_corr_stack_%s.mat', conditionID);
save(space_corr_stack_filename, 'corr_space');

mean_time_corr_filename=sprintf('mean_time_corr_%s.mat', conditionID);
save(mean_time_corr_filename, 'mean_corr_time');

mean_space_corr_filename=sprintf('mean_space_corr_%s.mat', conditionID);
save(mean_space_corr_filename, 'mean_corr_space');
% 
region_filename=sprintf('region_corr_avg_%s.csv', conditionID);
% save(region_filename, 'region_avg');

csvwrite(region_filename, region_avg);

end