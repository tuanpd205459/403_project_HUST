% draft main.m
%
clc; clear; close all;
load("main1.mat");
load("main_tai_tao_pha_bu.mat");
load('cat_anh.mat');

% 1. Tính toán sự khác biệt
result = reconSurface - Z_interp;

width_mm = 4.899;    % chiều rộng thực tế (mm)
height_mm = 3.660;   % chiều cao thực tế (mm)

% 2. Lấy kích thước ảnh
[rows, cols] = size(result);

% 3. Tạo lưới tọa độ
% Kích thước ảnh mong muốn
target_cols = 1061;
target_rows = 1421;

% Tạo lưới tọa độ mới (chuẩn hóa về kích thước mong muốn)
[X, Y] = meshgrid(...
    linspace(0, width_mm, cols), ...
    linspace(0, height_mm, rows));
%%
% 4. Vẽ ảnh 3D
figure;
surf(X, Y, reconSurface);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Nanomet');
title('Ảnh bề mặt sau tái tạo');
colorbar;
%%
figure;
surf(X, Y, Z_interp);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Nanomet');
title('Bề mặt nội suy');
colorbar;
%%
figure();
surf(X, Y, result);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Nanomet');
title('Bề mặt ảnh sau bù pha');
colorbar;
%%
save('bu_pha.mat');
