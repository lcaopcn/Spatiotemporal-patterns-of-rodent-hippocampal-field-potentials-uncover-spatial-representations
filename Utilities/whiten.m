function [X,wh,dw,zp] = whiten(X,fudgefactor,useEig)
if nargin < 3
    useEig = 1;
end
if nargin == 1
    fudgefactor = 0;
end
X = bsxfun(@minus, X, mean(X,2));
A = X*X'/size(X,2);
if useEig
[V,D] = eig(A);
%D = diag(max(diag(D),eps));
D = D + (fudgefactor-min(0,min(diag(D))-eps))*eye(size(D,1));
d = diag(D); %ones(size(D,1),1);%+1e-10
else
    [V,D] = svd(A,'econ');
    d = diag(D);
end
dsqrtinv = real(d + ones(size(d))*eps).^(-0.5);
wh = diag(dsqrtinv)*V';
dw = V*sqrt(D);%Ex * sqrt (Dx);
%zp = V*V';%*diag([ones(1,size(D,1)-1) 0])
zp = V*diag(dsqrtinv)*V';
%zerophaseMatrix = E*inv (sqrt (D))*E';
%zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';

X = wh*X;


% if ndims(im) == 3
%     im = rgb2gray(im);
% end
% 
% [n1 n2] = size(im);
% 
% [fx fy]=meshgrid(-n2/2:n2/2-1,-n1/2:n1/2-1);
% rho=sqrt(fx.*fx+fy.*fy);
% %f_0=0.4*N;
% %filt=rho.*exp(-(rho/f_0).^4);
% filt = rho.*exp(-.5*(rho/(.7*min(n1/2,n2/2))).^2);
% 
% If=fft2(im);
% imw=real(ifft2(If.*fftshift(filt)));
% %imwa=real(ifft2(If.*fftshift(rho)));
% 
% %imw = sqrt(0.1)*imw/sqrt(var(imw(:)));
% imw = imw/sqrt(var(imw(:)))*sqrt(var(im(:)));