function surface = reconstruct_wavefront(zernike_coeffs, order, grid_size, x_axis, y_axis)
    % Hàm tái tạo bề mặt sóng từ các hệ số đa thức Zernike
    % zernike_coeffs: Hệ số đa thức Zernike
    % order: Bậc của các đa thức Zernike
    % grid_size: Kích thước của lưới điểm trên mặt phẳng 
    % x_axis: Vector [xmin, xmax] xác định phạm vi trục x
    % y_axis: Vector [ymin, ymax] xác định phạm vi trục y

    % Kiểm tra đầu vào
    if length(x_axis) ~= 2 || length(y_axis) ~= 2
        error('x_axis và y_axis phải là vector có 2 phần tử [min, max].');
    end

    % Khởi tạo mặt phẳng lưới
    [X, Y] = meshgrid(linspace(x_axis(1), x_axis(2), grid_size), ...
                      linspace(y_axis(1), y_axis(2), grid_size));
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X);

    % Khởi tạo bề mặt sóng
    surface = zeros(size(X));

    % Lặp qua các bậc và hệ số để tái tạo bề mặt sóng
    index = 1;
    for n = 0:order
        for m = -n:2:n
            if index <= length(zernike_coeffs)
                % Lấy hệ số tương ứng
                Z = zernike_coeffs(index);
                
                % Tính giá trị Zernike polynomial
                Zmn = reconstructZernike.ZernikePolynomial(n, m, R, Theta);
                
                % Cộng dồn vào bề mặt sóng
                surface = surface + Z * Zmn;

                index = index + 1;
            end
        end
    end
end