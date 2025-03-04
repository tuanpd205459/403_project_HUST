
function Zmn = ZernikePolynomial(n, m, R, Theta)
    % Hàm tính Zernike polynomial
    % n: Bậc n của Zernike
    % m: Bậc m của Zernike
    % R: Bán kính
    % Theta: Góc theta
    
    % Tính radial Zernike polynomial
    RadialPoly = zeros(size(R));
    for k = 0:floor((n - abs(m)) / 2)
        c = (-1)^k * factorial(n - k) / ...
            (factorial(k) * factorial((n + abs(m)) / 2 - k) * factorial((n - abs(m)) / 2 - k));
        RadialPoly = RadialPoly + c * R.^(n - 2 * k);
    end
    
    % Tính Zernike polynomial
    if m >= 0
        Zmn = RadialPoly .* cos(m * Theta);
    else
        Zmn = RadialPoly .* sin(abs(m) * Theta);
    end
end
