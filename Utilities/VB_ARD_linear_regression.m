function [w, V, invV, logdetV, an, bn, E_a, L] =   VB_ARD_linear_regression(X, y)

% estimates w for a linear regression  y = Xw
% using sparse Bayesian regularization with automatic relevance determination.

% y and w are vectors, X are matrices

% The underlying generative model assumes p(y | x, w, tau) = N(y | w'x, tau^-1),
% with x and y being the rows of the given X and y. w and tau are
% assigned the conjugate normal inverse-gamma prior 
% p(w, tau | alpha) =N(w | 0, (tau A)^-1) Gam(tau | a0, b0), with A being a diagonal matrix
% with the vector alpha = (alpha_1, ..., alpha_D)' along its diagonal,
% and the hyperpriors   p(alpha_i) = p(alpha | c0, d0). 

% The prior parameters a0, b0, c0, d0 are set such that the prior is non-informative.
% The returned posterior (computed by variational Bayesian inference) is of
% the form N(w1 | w, tau^-1 V) Gam(tau | an, bn). 

% Also, the mean vector E_a = E(alpha) is returned, together with the inverse of V, and its log-determinant. 
% L is the variational bound of the model, and is a lower bound on the log-model evidence ln p(y | X).

% assign uninformative priors
a0 = 1e-2;
b0 = 1e-4;
c0 = 1e-2;
d0 = 1e-4;

[N, D] = size(X);
% pre-process data
X_corr = X' * X;
Xy_corr = X' * y;
an = a0 + N / 2;
cn = c0 + 1 / 2;

% iterate to find hyperparameters
L_last = -realmax;
max_iter = 1000;
E_a = ones(D, 1) * c0 / d0;

for iter = 1:max_iter
    
    % covariance and weight of linear model
    invV = diag(E_a) + X_corr;
    V = inv(invV);
    logdetV = - logdet(invV);
    w = V * Xy_corr;

    % parameters of noise model (an remains constant)
    sse = sum((X * w - y) .^ 2);
    bn = b0 + 0.5 * (sse + sum(w .^ 2 .* E_a));
    E_t = an / bn;

    % hyperparameters of covariance prior (cn remains constant)
    dn = d0 + 0.5 * (E_t .* w .^ 2 + diag(V));
    E_a = cn ./ dn;

    % variational bound
    L = - 0.5 * (E_t * sse + sum(sum(X .* (X * V)))) + 0.5 * logdetV ...
          - b0 * E_t + gammaln(an) - an * log(bn) + an + D * gammaln(cn) - cn * sum(log(dn));

    % variational bound must grow!
    if L_last > L
        fprintf('Last bound %6.6f, current bound %6.6f\n', L_last, L);
        error('Variational bound should NOT reduce');
    end

    % stop if change in variation bound is < 0.001%
    if abs(L_last - L) < abs(0.00001 * L)
        break
    end
    L_last = L;
end

if iter == max_iter
     warning('VB: maxIter', ' reached maximum number of iterations.');
end

% augment variational bound with constant terms
L = L - 0.5 * (N * log(2 * pi) - D) - gammaln(a0) + a0 * log(b0) + D * (- gammaln(c0) + c0 * log(d0));
