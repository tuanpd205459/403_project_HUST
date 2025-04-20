function surface = zernike_reconstruction(coefficients, real_radius, lambda, surface_title)
    % Hàm tái tạo bề mặt từ hệ số Zernike
    % coefficients: vector hệ số Zernike
    % real_radius: bán kính khẩu độ (mm)
    % lambda: bước sóng ánh sáng (nm)
    % surface_title: tiêu đề của hình ảnh

    num_coeff = length(coefficients);
    grid_size = 1061;

    [x, y] = meshgrid(linspace(-real_radius, real_radius, grid_size));
    [theta, rho] = cart2pol(x, y);
    rho = rho / real_radius;   % Chuẩn hóa về [0,1]
    mask = rho <= 1;

    if nargin < 4 || isempty(surface_title)
        surface_title = "Bề mặt tái tạo từ đa thức Zernike";
    end

    surface = zeros(size(x));
    n = 0; m = 0;

    for k = 1:num_coeff
        Z = zernike_polynomial(n, m, rho(mask), theta(mask));
        surface(mask) = surface(mask) + coefficients(k) * Z;

        if m == n
            n = n + 1;
            m = -n;
        else
            m = m + 2;
        end
    end

    surface = surface * lambda;       % Đơn vị: nanomet
    surface(~mask) = NaN;

    % Hiển thị bề mặt 3D
    figure;
    surf(x, y, surface, 'EdgeColor', 'none');
    title(surface_title, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Độ cao sóng (nm)');
    colormap(jet);
    colorbar;
    view(3);
%     axis equal;
%     shading interp;
%     camlight headlight;
%     lighting phong;
end


%% Hàm Zernike cơ bản
function Z = zernike_polynomial(n, m, rho, theta)
    R = zernike_radial(n, abs(m), rho);
    if m >= 0
        Z = R .* cos(m * theta);
    else
        Z = R .* sin(abs(m) * theta);
    end
end

function R = zernike_radial(n, m, rho)
    R = zeros(size(rho));
    for s = 0:((n - m)/2)
        c = (-1)^s * factorial(n - s) / ...
            (factorial(s) * factorial((n + m)/2 - s) * factorial((n - m)/2 - s));
        R = R + c * rho.^(n - 2*s);
    end
end
