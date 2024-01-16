clc
clear all
% rng('default')

% tol = 1e-3;
tol = {1e-3, 1e-4};

n = [100, 100, 100];
r = [5, 5, 5];
d = size(n,2);

%%% Tensor Data
C = tenrandn(r);
U = cell(d,1);
for i=1:d
    [U{i}, ~] = qr(randn(n(i),r(i)), 0);
end
X = ttm(C,U);
clear C U

Noise = tenrandn(n);
X = X/norm(X) + .1*Noise/norm(Noise);
clear Noise

%%% Initial U
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(n(i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~, ~] = RPCD(X,r,'tol',tol,'init',Uinit);

% RPCD plus
[~, ~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

