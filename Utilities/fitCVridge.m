% Fit a regularized GLM (input X, output y) using n-fold Cross Validation
% Returns a structure m, with params and log-likelihood information
% m has fields...
%   b       [p x nfoldcv x length(nu_vec)] matrix of parameters
%   breg    [p x nfoldcv]                  MAP estimates at the optimized hyperparam
%   llhdr       training-set negative log-likelihoods (per bin)
%   llhdt       test-set negative log-likelihoods (per bin)
%   llhdt0      test-set negative log-likelihoods (per bin) for homogeneous process
%   llrt        test-set log likelihood ratios relative to homogeneous process
%   numspks     test-set spikes to get llrt/spk
%   numbins     test-set bins to get llrt/s
%   exitflag    output from fminunc to determine convergence

function m = fitCVridge(X,y,nfoldcv,nu_vec)

% Don't penalize the baseline parameter...
penalty = ones(size(X,2)+1,1);
penalty(1) = 0;

% Add intercept...
X = [X(:,1)*0+1 X];

% Cross-Validation
mcvpcv = getCVidx(size(y,1),10);
for i=1:nfoldcv
    fprintf('iF %02i/%02i,',i,nfoldcv)
    if nfoldcv==1
        idx_tr = 1:size(y,1);
        idx_ts = 1:size(y,1);
    else
        idx_tr = mcvpcv.tr{i};
        idx_ts = mcvpcv.ts{i};
    end
    idx_tr = idx_tr(isfinite(y(idx_tr)));
    idx_ts = idx_ts(isfinite(y(idx_ts)));
    
    m.cvidx_tr{i} = idx_tr;
    m.cvidx_ts{i} = idx_ts;
    
    % Fit a homogeneous model to compare others against...
    bmu = mean(y(idx_tr));
    m.llhdt0(i) = sum((y(idx_ts)-bmu).^2)/length(idx_ts);
    m.bmu(i) = bmu;
    
    % Evaluate the saturated model...
    m.llhdts(i) = 0;
    C = X(idx_tr,:)'*X(idx_tr,:);
    
    % Loop over regularization hyper-parameters...
    for j=1:length(nu_vec)
        m.b(:,i,j)  = inv(C+nu_vec(j)*diag(penalty))*X(idx_tr,:)'*y(idx_tr);
        m.llhdr(i,j) = sum((y(idx_tr)-X(idx_tr,:)*m.b(:,i,j)).^2)/length(idx_tr);
        m.llhdt(i,j) = sum((y(idx_ts)-X(idx_ts,:)*m.b(:,i,j)).^2)/length(idx_ts);
        m.llrt(i,j)  = (-m.llhdt(i,j)+m.llhdt0(i))/log(2);
        m.pseudoR2(i,j) = 1-(m.llhdts(i)-m.llhdt(i,j))/(m.llhdts(i)-m.llhdt0(i));
        
        if (penalty'*m.b(:,i,j))==0 
            fprintf('Zero reached...');
            break;
        end
    end
    if isfield(m,'llrt')
        [m.llrtmax(i),id]=max(m.llrt(i,:));
        m.breg(:,i) = m.b(:,i,id);
    end
    
    % Estimated (Cross-Validated) Firing Rate...
    m.lambda(idx_ts) = X(idx_ts,:)*m.breg(:,i);
    % (Cross-Validated) Cross-Correlation...
    m.xcorr(i) = diag(corrcoef(m.lambda(idx_ts),y(idx_ts)),1);
    m.xcorr_rect(i) = diag(corrcoef(max(0,m.lambda(idx_ts)),y(idx_ts)),1);
end