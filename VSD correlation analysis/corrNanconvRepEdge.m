function Dout = corrNanconvRepEdge(D,kernel)
% nanconvn - performs n-dimensional convolution, ignoring NaNs
% and weighting each point by the actual number of datapoints involved.
% zero-padding 


Dpad = padarray(D,round(size(kernel)/2)-1,'replicate','both'); % length of convolution is MAX([LENGTH(A)+LENGTH(B)-1,LENGTH(A),LENGTH(B)].  in our case this will almost certainly be length(a)+length(b)-1

%create a mask of non-NaN values
mask = ~isnan(Dpad);
%convolving the mask with B will give weights representing how many and how much
%actual data points actually contributed to the convolved values 

%replace NaNs by zeros in A and B
Dpad(isnan(Dpad))=0;

%convolve the NaN-less version of A with B and divide by the weights
Dout = convn(Dpad,kernel,'valid')./convn(mask,kernel,'valid')*sum(kernel(:)); % sum(kernel(:)) gives the sum of all elements in kernel - this makes the output have the same value as output from convRepEdge, for non-NaN inputs

%Note that some points will be undefined (those that don't receive any contribution from any point of A)

