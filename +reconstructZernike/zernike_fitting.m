function zernike_fitting(reconsurface)
    % Load surface data from the input variable
    Z = reconsurface;
    [rows, cols] = size(Z);

    % Generate X and Y coordinates
    [X, Y] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));

    % Convert Cartesian coordinates to polar coordinates
    [theta, rho] = cart2pol(X, Y);

    % Maximum order of Zernike polynomials
    maxOrder = 10;

    % Compute Zernike polynomials and fit the surface
    ZernikeCoeffs = fitZernike(rho, theta, Z, maxOrder);

    % Reconstruct the surface using Zernike polynomials
    Z_fit = reconstructZernike(rho, theta, ZernikeCoeffs);

%     % Plot original surface
%     figure;
%     surf(X, Y, Z);
%     title('Original Surface');
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');

    % Plot fitted surface
    figure;
    surf(X, Y, Z_fit);
    title('Fitted Surface with Zernike Polynomials');
    xlabel('X');
    ylabel('Y');
    zlabel('Z_ fit');
shading interp;
colormap('jet');

colorbar;
end

% Function to fit Zernike polynomials to the surface
function ZernikeCoeffs = fitZernike(rho, theta, Z, maxOrder)
    numTerms = (maxOrder + 1) * (maxOrder + 2) / 2;
    ZernikeCoeffs = zeros(numTerms, 1);

    % Generate Zernike polynomials and solve for coefficients
    index = 1;
    for n = 0:maxOrder
        for m = -n:2:n
            ZernikePoly = zernike(n, m, rho, theta);
            ZernikeCoeffs(index) = sum(sum(Z .* ZernikePoly)) / sum(sum(ZernikePoly .* ZernikePoly));
            index = index + 1;
        end
    end
end

% Function to reconstruct surface using Zernike polynomials
function Z_fit = reconstructZernike(rho, theta, ZernikeCoeffs)
    maxOrder = floor(sqrt(2 * length(ZernikeCoeffs)) - 1);
    Z_fit = zeros(size(rho));

    index = 1;
    for n = 0:maxOrder
        for m = -n:2:n
            ZernikePoly = zernike(n, m, rho, theta);
            Z_fit = Z_fit + ZernikeCoeffs(index) * ZernikePoly;
            index = index + 1;
        end
    end
end

% Function to compute Zernike polynomials
function Z = zernike(n, m, rho, theta)
    if m > 0
        Z = sqrt(2) * radialZernike(n, m, rho) .* cos(m * theta);
    elseif m < 0
        Z = sqrt(2) * radialZernike(n, -m, rho) .* sin(-m * theta);
    else
        Z = radialZernike(n, 0, rho);
    end
end

% Function to compute radial Zernike polynomials
function R = radialZernike(n, m, rho)
    R = zeros(size(rho));
    for s = 0:((n - abs(m)) / 2)
        c = ((-1)^s * factorial(n - s)) / (factorial(s) * factorial((n + abs(m)) / 2 - s) * factorial((n - abs(m)) / 2 - s));
        R = R + c * rho.^(n - 2 * s);
    end
end