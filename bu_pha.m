% draft main.m
%
clc; clear; close all;
load("reconSurface.mat");
load("averageMatrix.mat");
load('data_saunoisuy.mat');

% 1. Tính toán sự khác biệt
result = reconSurface - Z_interp;

% 2. Lấy kích thước ảnh
[rows, cols] = size(result);

% 3. Tạo lưới tọa độ
[X, Y] = meshgrid(1:cols, 1:rows);
%%
% 4. Vẽ ảnh 3D
figure;
surf(X, Y, reconSurface);
shading interp;
colormap('jet');
xlabel('X (pixel)');
ylabel('Y (pixel)');
zlabel('Cường độ');
title('ban dau');
colorbar;
%%
% 4. Vẽ ảnh 3D
figure;
surf(X, Y, Z_interp);
shading interp;
colormap('jet');
xlabel('X (pixel)');
ylabel('Y (pixel)');
zlabel('Cường độ');
title('Pha bu');
colorbar;
%%
% 4. Vẽ ảnh 3D
figure;
surf(X, Y, result);
shading interp;
colormap('jet');
xlabel('X (pixel)');
ylabel('Y (pixel)');
zlabel('Cường độ');
title('Bề mặt ảnh sau bu pha');
colorbar;
save('bu_pha.mat');
