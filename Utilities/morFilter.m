function Xf = morFilter(X,Fc,Fb,Fs,norm,win)

%X = t x n matrix, where n = # channels, t = # timesteps
%Fc = center frequency
%Fs = sampling rate
%norm = normalize output of each channel to 1 (default 0)
%win = window for filter (default 1 s)
%Xf = complex-valued filtered output
%Fb - This parameter controls bandwidth, may need to be changed for your purposes
%
%Fc = 8.4 for hippocampus
%Fb = 1/500, win = 1 for hippocampal LFP
%Fb = 1/5000, Fc = 160, win = 5 for ripples


if ~exist('norm','var')
    norm = 0;
end
if ~exist('win','var')
    win = 1;
end

[psi,x] = cmorwavf(-win,win,Fs*2*win,Fb,Fc);

Xf = zeros(size(X));
for j = 1:size(X,2)
    Xf(:,j) = flipud(filter(conj(psi),1,flipud(filter(psi,1,X(:,j)))));
    if norm, Xf(:,j) = Xf(:,j)/std(Xf(:,j)); end
end