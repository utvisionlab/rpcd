% load_boats:  load coil-100 dataset
%
% This is a code for "D-TUCKER: Fast and Memory-Efficient Tucker Decomposition for Dense Tensors", submitted to ICDE 2020.
% Authors: Jun-Gi Jang (elnino9158@gmail.com), Seoul National University
%          U Kang (ukang@snu.ac.kr), Seoul National University
%
% Return
% X: [128 128 72 100] tensor
%

function X = load_coil()

X = imread('./data/coil-100/obj1__0.png');
X = rgb2gray(X);
for i=2:72
   tmp_path = strcat('./data/coil-100/obj1__', ...
                    int2str((i-1)*5), '.png');
   tmp = imread(tmp_path);
   tmp = rgb2gray(tmp);
   X = cat(3, X, tmp);
end

path = './data/coil-100/obj';
for i=2:100
    Y = imread(strcat(path, int2str(i), '__0.png'));
    Y = rgb2gray(Y);
    for j=2:72
        tmp_path = strcat(path, int2str(i), '__', int2str((j-1)*5), '.png');
        tmp = imread(tmp_path);
        tmp = rgb2gray(tmp);
        Y = cat(3, Y, tmp);
    end
    X = cat(4, X, Y);
end
end
