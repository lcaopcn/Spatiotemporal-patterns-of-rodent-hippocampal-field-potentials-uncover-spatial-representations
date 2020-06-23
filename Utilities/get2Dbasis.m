function [b,Xrbf] = get2Dbasis(btype,n,X,s)

btype=lower(btype);
b.n=n;
b.type=btype;
if strcmp(btype,'gaussian')
    b.mxv = linspace(0,1,n(1));
    b.myv = linspace(0,1,n(2));
    [xm,ym] = meshgrid(b.mxv,b.myv);
    b.mxy = [xm(:) ym(:)];
    if nargin>3
        b.s=s;
    else
        b.s = diag([mean(diff(b.mxv)) mean(diff(b.myv))]).^2;
    end
    
    [xm,ym]=meshgrid(linspace(0,1,256),linspace(0,1,256));
    b.x0=[xm(:) ym(:)];
    b.rbf_basis = get2Dbasis_gaussian(b.x0,b.mxy,b.s);
    if nargin>2
        Xrbf = get2Dbasis_gaussian(X,b.mxy,b.s);
    end
end

c=1;
for i=1:n(1)
    for j=1:n(2)
        subplot_tight(n(1),n(2),c)
        imagesc(reshape(b.rbf_basis(:,c),[256 256]))
        axis off
        axis image
        c=c+1;
    end
end


function Xrbf = get2Dbasis_gaussian(X,rbfm,rbfs)

Xrbf=[];
for i=1:size(rbfm,1)
    Xrbf(:,i) = exp(-sum(bsxfun(@times,bsxfun(@minus,X,rbfm(i,:))*inv(rbfs),bsxfun(@minus,X,rbfm(i,:))),2)/2)/sqrt(det(rbfs))/(2*pi);
end