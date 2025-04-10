% Hàm nội suy bề mặt 2D bằng đa thức Zernike trong MATLAB
%-------------------------------------------------------------
function zernike_fitting2(surface, resolution, max_order)
    [rows, cols] = size(surface);
    [X, Y] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
    coefficients = zeros((max_order + 1) * (max_order + 2) / 2, 1);

    % Tính toán hệ số Zernike
    index = 1;
    for n = 0:max_order
        for m = -n:2:n
            Z = zernike_polynomial(n, m, X, Y);
            C = sum(surface .* Z, 'all') / sum(Z.^2, 'all');
            coefficients(index) = C;
            index = index + 1;
        end
    end

    % Nội suy bề mặt bằng đa thức Zernike
    [X, Y] = meshgrid(linspace(-1, 1, resolution), linspace(-1, 1, resolution));
    Z = zeros(size(X));
    index = 1;
    for n = 0:max_order
        for m = -n:2:n
            Z = Z + coefficients(index) * zernike_polynomial(n, m, X, Y);
            index = index + 1;
        end
    end

    % Cắt thành hình tròn
    Z = mask_circle(Z);

    % Hiển thị kết quả
    figure;
    surf(X, Y, Z);
    title('Nội suy bề mặt 2D bằng đa thức Zernike - Hình tròn');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;

end

% Hàm tính giá trị đa thức Zernike
function Z = zernike_polynomial(n, m, x, y)
    if mod(n-m, 2) ~= 0 || abs(m) > n
        Z = zeros(size(x));
        return;
    end

    [theta, rho] = cart2pol(x, y);
    rho = min(rho, 1);

    R = 0;
    for s = 0:((n-abs(m))/2)
        c = (-1)^s * factorial(n-s) / (factorial(s) * factorial((n+abs(m))/2-s) * factorial((n-abs(m))/2-s));
        R = R + c * rho.^(n-2*s);
    end

    if m >= 0
        Z = R .* cos(m * theta);
    else
        Z = R .* sin(abs(m) * theta);
    end
end

% Hàm cắt thành hình tròn
function Z = mask_circle(Z)
    [rows, cols] = size(Z);
    [X, Y] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
    mask = sqrt(X.^2 + Y.^2) <= 1;
    Z(~mask) = NaN;
end


