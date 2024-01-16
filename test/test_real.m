clc
clear all

% tol = 1e-3;
tol = {1e-3, 1e-4};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Brainq data 
fprintf('\nRead brainq data\n');
path = './data/brainq/BRAINQ.tensor';
X = import_data(path);
X = tensor(X);
fprintf('Reading brainq data is done\n');
r = [10 10 5];

d = ndims(X);
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(size(X,i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~,~] = RPCD(X,r,'tol',tol,'init',Uinit);
% RPCD plus
[~,~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Air quality data
fprintf('\nRead Air quality data\n');
path = './data/airquality/airquality.tsv';
X = import_data(path);
X = tensor(X);
fprintf('Reading air quality data is done\n');
r = [10 10 5];

d = ndims(X);
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(size(X,i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~,~] = RPCD(X,r,'tol',tol,'init',Uinit);
% RPCD plus
[~,~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HSI data
fprintf('\nRead HSI data\n');
X = cell([1 8]);
dims_tmp = zeros([8 3]);
for i=1:8
    path = './data/hsi/input/data';
    path = strcat(path, int2str(i), '.mat');
    tmp = load(path);
    X{1,i} =cell2mat(struct2cell(tmp));
    dims_tmp(i,:) = size(X{i});
end
maxlength = max(dims_tmp, [], 1);
for i = 1:8
    tmp_dims = size(X{i});
    for j =1:3
        curr_dims = tmp_dims(j);
        tmp_dims(j) = maxlength(j);
        add_dims = tmp_dims;
        add_dims(j) = tmp_dims(j)-curr_dims;
        X{i} = cat(j, X{i}, zeros(add_dims));
    end
end
X = cat(4, X{1}, X{2}, X{3}, X{4}, X{5}, X{6}, X{7}, X{8});
X = tensor(X);
fprintf('Reading HSI data is done\n');
r = [10 10 10 5];

d = ndims(X);
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(size(X,i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~,~] = RPCD(X,r,'tol',tol,'init',Uinit);
% RPCD plus
[~,~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coil data
fprintf('\nRead coil-100 data\n');
X = tensor(double(load_coil()));
fprintf('Reading coil-100 data is done\n');
r = [5 5 5 5];

d = ndims(X);
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(size(X,i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~,~] = RPCD(X,r,'tol',tol,'init',Uinit);
% RPCD plus
[~,~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Yale face data
fprintf('\nRead Yale face data\n');
load("Yale_64x64.mat")
data = fea(:,:);
X = tenzeros([64,64,11,15]);
for i=1:15
    for j=1:11
        X(:,:,j,i) = reshape(data((i-1)*11+j,:), [64,64]);
    end
end
fprintf('Reading Yale face data is done\n');
r = [5 5 5 5];

d = ndims(X);
Uinit = cell(d,1);
for i=1:d
%     [Uinit{i}, ~] = qr(randn(size(X,i),r(i)), 0);
    Uinit{i} = eye(size(X,i),r(i));
end

% RPCD
[~,~] = RPCD(X,r,'tol',tol,'init',Uinit);
% RPCD plus
[~,~] = RPCD_plus(X,r,'tol',tol,'init',Uinit);

