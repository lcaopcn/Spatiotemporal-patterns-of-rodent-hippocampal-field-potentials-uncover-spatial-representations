
%% Plot Decoding Results...
idx=1:size(f,1);
[tmp,maxpost]=max(f');
xyhat = pvec(maxpost,:);
idx=1:size(xyhat,1);

figure(1)
clf
subplot(2,1,1)
plot(rpos(:,1),'k.')
hold on
plot(xyhat(:,1),'r.')
hold off
axis tight
ylabel('X-Position')

subplot(2,1,2)
plot(rpos(:,2),'k.')
hold on
plot(xyhat(:,2),'r.')
hold off
axis tight
ylabel('Y-Position')
xlabel('Time [bins]')

%% Absolute error...

err = sqrt(sum(((rpos(idx,:)-xyhat).*[xlen,ylen]).^2,2));

figure(3)
hist(abs(err),50);
axis tight
xlabel('Decoding Error [cm]')
ylabel('Frequency')