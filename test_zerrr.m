% Thiết lập lưới điểm
grid_size = 100;
x = linspace(-1, 1, grid_size);
y = linspace(1, -1, grid_size);
[X, Y] = meshgrid(x, y);

% Tạo hàm mẫu z(x,y)
z = X.^2 + Y.^3;
figure;
subplot(2, 2, 1);
imagesc(z);
colormap jet;
colorbar;
title('Hàm mẫu z(x,y)');

% Thiết lập danh sách các chỉ số (n, m) của đa thức Zernike
indices = [];
for n = 0:3
    for m = -n:2:n
        indices = [indices; n m];
    end
end

% Tính hệ số Zernike
moments = zernike_moments(z, indices);

% Hiển thị hệ số Zernike
disp('    Coeff     n        m');
disp([moments indices]);

% Hiển thị các đa thức Zernike
zernike_sum = zeros(grid_size, grid_size);
zern_mats = zernike_mats(zernike_sum, indices);

% Vẽ các đa thức Zernike
for i = 1:min(3, size(zern_mats, 3))
    subplot(2, 2, i + 1);
    imagesc(zern_mats(:, :, i));
    colormap jet;
    colorbar;
    title(['Đa thức Zernike n=' num2str(indices(i, 1)) ', m=' num2str(indices(i, 2))]);
end
