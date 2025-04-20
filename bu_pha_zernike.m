% draft main.m
%
clc; clear;

load("main1.mat");
load("main_tai_tao_pha_bu.mat");
load('cat_anh.mat');
load('main1_zernike.mat');

%%


% 1. Tính toán sự khác biệt
reconSurface = reconSurface(1:1061, :);     % cắt ảnh hình tròn 1061 x1061;
result = reconSurface - surface_zernike3;

% width_mm = 4.899;    % chiều rộng thực tế (mm)
height_mm = 3.660;   % chiều cao thực tế (mm)
width_mm = 3.660;

% 2. Lấy kích thước ảnh
[rows, cols] = size(result);

% 3. Tạo lưới tọa độ
% Kích thước ảnh mong muốn
target_cols = 1061;
target_rows = 1061;

% Tạo lưới tọa độ mới (chuẩn hóa về kích thước mong muốn)
[X, Y] = meshgrid(...
    linspace(0, width_mm, cols), ...
    linspace(0, height_mm, rows));
%%
[reconSurface, dimensional_reconSurface] = myConvertUnit(reconSurface);
[result, dimensional_result] = myConvertUnit(result);
[surface_zernike3, dimensional_surface_zernike3] = myConvertUnit(surface_zernike3);


% 4. Vẽ ảnh 3D
figure;
surf(X, Y, reconSurface);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel(dimensional_reconSurface);
title('Ảnh bề mặt sau tái tạo (Zernike)');
colorbar;
%%
figure;
surf(X, Y, surface_zernike3);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel(dimensional_surface_zernike3);
title('Bề mặt tai tao bang Zernike');
colorbar;
%%
figure();
surf(X, Y, result);
shading interp;
colormap('jet');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Nanomet');
title('Bề mặt ảnh sau bù Zernike');
zlabel(dimensional_result);
colorbar;
%%
save('bu_pha_zernike.mat');


function [reconSurface, dimensional] = myConvertUnit(reconSurface)
    % Kiểm tra giá trị lớn nhất để xác định đơn vị
%     if (max(reconSurface(:))-min(reconSurface(:))) < 10
    if max(abs(reconSurface(:))) > 1000
        scaleFactor = 10^3; % Chuyển từ nanomet sang micromet
        dimensional = 'micromet';
    else
        scaleFactor = 1;  
        dimensional = 'nanomet';
    end

    % Chuyển đổi đơn vị
    reconSurface = reconSurface / scaleFactor;
end


