
function Z_flat =  visualize_surface_flat_once(Z)
% Loại nghiêng chính xác bằng cách chiếu vuông góc 1 lần duy nhất lên mặt phẳng chuẩn

% 1. Tạo bề mặt rõ nghiêng và gồ ghề

Z = imgaussfilt(Z, 1);  % làm mịn nhẹ
[x,y ] = size(Z);

[x, y] = meshgrid(1:x, 1:y);


% 2. Fit mặt phẳng: z = ax + by + c
A = [x(:), y(:), ones(numel(x),1)];
coeffs = A \ Z(:);
a = coeffs(1); b = coeffs(2); c = coeffs(3);

% 3. Tính độ lệch vuông góc đến mặt phẳng
denom = sqrt(a^2 + b^2 + 1);
Z_flat = (Z - (a*x + b*y + c)) / denom;

% 4. Vẽ bề mặt ban đầu và sau hiệu chỉnh
figure;
subplot(1,2,1); surf(x, y, Z, 'EdgeColor', 'none');
title('Bề mặt ban đầu (nghiêng + gồ)'); xlabel('x'); ylabel('y'); zlabel('Z');
axis tight; view(3); colorbar;

subplot(1,2,2); surf(x, y, Z_flat, 'EdgeColor', 'none');
title('Bề mặt sau khi loại nghiêng (chiếu 1 lần)'); xlabel('x'); ylabel('y'); zlabel('Z flat');
axis tight; view(3); colorbar;

% 5. Vẽ mặt cắt qua y = giữa
idc = round(size(y,1)/2);
figure;
plot(x(idc,:), Z(idc,:), 'b', 'LineWidth', 1.5); hold on;
plot(x(idc,:), Z_flat(idc,:), 'r--', 'LineWidth', 1.5);
legend('Gốc', 'Sau loại nghiêng');
title('Mặt cắt ngang qua y = 50'); xlabel('x'); ylabel('Z'); grid on;
end
