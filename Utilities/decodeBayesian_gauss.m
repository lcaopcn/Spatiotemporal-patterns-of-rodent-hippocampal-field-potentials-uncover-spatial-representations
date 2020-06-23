
function pall = decodeBayesian_gauss(y,lam,sigma,od)

if nargin<4, od=0; end

pall=zeros(size(y,1),size(lam,1));
for i=1:size(y,1)
    for p=1:size(lam,1)
        P(:,p) = normpdf(y(i,:)',lam(p,:)',sigma)+od;
    end
    P = bsxfun(@rdivide,P,nansum(P,2));
    pall(i,:) = nansum(log(P));
end

pall=(exp(bsxfun(@minus,pall,max(pall))));
pall=bsxfun(@rdivide,pall,nansum(pall,2));

pall(~isfinite(nansum(pall,2)),:)=1/size(pall,2);
pall(isnan(pall))=0;