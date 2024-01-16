function [T, Error, time_stamp] = RPCD(X,r,varargin)
% Riemannian Preconditioned Coordinate Descent
%
% T = RPCD(X,r) compute the rank-(r1,r2,...,rd) approximation of the
% dense tensor X according to the specified dimensions in vector r. 
% The result is the ttensor T and the relative error Error.
%
% T = RPCD(X,r,'param',value,...) specifies optional parameters and values.
% Valid parameters and their default values are:
%       'tol'      :  Tolerence in difference in relative error {1e-3}
%       'alpha'    :  Step-size length {1}
%       'maxiters' :  Maximum number of iterations {20}
%       'init'     :  Initial guess {cell array}
%

% Extract the dimension informations of tensor X
n = size(X);
d = ndims(X);

% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('tol',1e-3,@(x) (isscalar(x) || iscell(x)));
params.addParameter('alpha',1,@isscalar);
params.addParameter('maxiters',20,@(x) isscalar(x) & x > 0);
params.addParameter('init', 'random', @(x) iscell(x));
params.parse(varargin{:});

% Copy from params object
tol = params.Results.tol;
alpha = params.Results.alpha;
maxiters = params.Results.maxiters;
init = params.Results.init;

if isscalar(tol)
    tol_outer = tol;
else
    tol_outer = tol{1};
end

% Initialization
if iscell(init)
    U = init;
else
    if strcmp(init,'random')
        fprintf('\nRANDOM INITIALIZATION');
        U = cell(d,1);
        for i = 1:d
            [U{i}, ~] = qr(randn(n(i),r(i)), 0);
        end
    end
end
C = ttm(X, U, 't');

% initial error
normX = norm(X);
normC = norm(C);
relErr = sqrt(normX^2 - normC^2)/normX;
time_stamp = 0;

fprintf('\nRPCD with Relative Error Delta:\n');

% Main loop
total_time = 0;
for k=1:maxiters
    tstart = tic;
    % Iterate over all d modes of the tensor
    for i=1:d
        Y = ttm(X, U, -i, 't');
        Yi = double(tenmat(Y, i));
        
        % Riemannian gradient
        G = - Yi*(Yi'*U{i});
        RG = U{i} + G;
        
        % Update
        dir = -RG;
        U{i} = U{i} + alpha*dir;
        % Retraction : Thin QR factorization
        [U{i}, ~] = qr(U{i}, 0); 
    end
    tElapsed = toc(tstart);
    total_time = total_time + tElapsed;
    time_stamp = [time_stamp; total_time];
    
    % Computing relative error
    normC = norm(U{i}'*Yi, 'fro');
    relErr = [relErr; sqrt(normX^2 - normC^2)/normX];
    
    % Stopping criteria
    relErrDelta = abs(relErr(end-1) - relErr(end));
    fprintf(' Iter %d: relErrDelta: %.2e relErr: %.5f in %.2f second\n', k, relErrDelta, relErr(end), tElapsed)
    
    % Error = relErr(end);
    Error = relErr;

    % Check for convergence
    if relErrDelta < tol_outer
        fprintf('Total time: %.2f second \n', total_time);
        break
    end
end

C = ttm(X, U, 't');
T = ttensor(C,U);

end

