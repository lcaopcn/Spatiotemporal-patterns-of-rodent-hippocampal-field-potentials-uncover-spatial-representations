function fs = smoothingFilter(f,S)

fs=f*0;
fs(1,:)=f(1,:);
for i=2:size(f,1)
    pri = fs(i-1,:)*S;
    pri = pri/nansum(pri);
    fs(i,:)=pri.*f(i,:);
end
fs = bsxfun(@rdivide,fs,nansum(fs,2));