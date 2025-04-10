% Script: cat_noisuy_anh.m
% Mục tiêu: Cắt ảnh thực 5.4x4.2 mm về 4.899x3.660 mm, làm mịn, nội suy lên 1421x1061
clc; clear; close all;
% Script: cat_noisuy_anh.m
% Mục tiêu: Cắt ảnh thực 5.4x4.2 mm về 4.899x3.660 mm, nội suy lên 1421x1061

%% 1. Đọc dữ liệu ảnh từ Excel
load('main_tai_tao_pha_bu.mat');
Z = averageMatrix;
[M, N] = size(Z);

% Kiểm tra dữ liệu hợp lệ
if isempty(Z) || M == 0 || N == 0
    error('Không đọc được dữ liệu từ file Excel');
end

%% 2. Kích thước pixel ban đầu
width_mm = 5.4;
height_mm = 4.2;
dx = width_mm / N;
dy = height_mm / M;

%% 3. Cắt ảnh về 4.899 x 3.660 mm
width_cut_mm = 4.899;
height_cut_mm = 3.660;
cols_cut = floor(width_cut_mm / dx);
rows_cut = floor(height_cut_mm / dy);

start_col = round((N - cols_cut) / 2) + 1;
start_row = round((M - rows_cut) / 2) + 1;

% Kiểm tra hợp lệ trước khi cắt
if start_col < 1 || start_row < 1 || ...
   start_col + cols_cut - 1 > N || ...
   start_row + rows_cut - 1 > M
    error('Chỉ số cắt ảnh không hợp lệ!');
end

% Cắt ảnh
Z_cut = Z(start_row:start_row + rows_cut - 1, start_col:start_col + cols_cut - 1);

%% 4. Nội suy lên lưới 1421 x 1061
target_rows = 1421;
target_cols = 1061;

[X_cut, Y_cut] = meshgrid(...
    linspace(0, width_cut_mm, cols_cut), ...
    linspace(0, height_cut_mm, rows_cut));

[X_new, Y_new] = meshgrid(...
    linspace(0, width_cut_mm, target_cols), ...
    linspace(0, height_cut_mm, target_rows));

Z_interp = interp2(X_cut, Y_cut, Z_cut, X_new, Y_new, 'linear');

%% 
% Chuyển sang đơn vị nanomet
Z_interp = Z_interp * 633;

%% 5. Vẽ ảnh 3D
figure;
surf(X_new, Y_new, Z_interp);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Nanomet');
title('Bề mặt ảnh sau cắt và nội suy');
colorbar;

%% 7. Lưu ảnh kết quả ra file Excel

save("cat_anh.mat");
