function zernike_reconstruction_non_circle(coefficients, width, height)
    % Function to reconstruct a surface from Zernike coefficients over a rectangular grid
    % coefficients: Vector of Zernike coefficients
    % width: Width of the rectangular grid
    % height: Height of the rectangular grid

    % Number of coefficients
    num_coeff = length(coefficients);

    % Generate a grid of points within the rectangle
    [x, y] = meshgrid(linspace(-width/2, width/2, 400), linspace(-height/2, height/2, 400));
    [theta, rho] = cart2pol(x, y);
    
    % Initialize the surface
    surface = zeros(size(x));

    % Calculate the surface using the Zernike coefficients
    n = 0; % Radial degree
    m = 0; % Azimuthal frequency
    for k = 1:num_coeff
        % Calculate the Zernike polynomial
        Z = zernike_polynomial(n, m, rho, theta);
        
        % Add the contribution of the current Zernike polynomial
        surface = surface + coefficients(k) * Z;

        % Update the radial degree and azimuthal frequency
        if m == n
            n = n + 1;
            m = -n;
        else
            m = m + 2;
        end
    end

    % Plot the reconstructed surface
    figure;
    surf(x, y, surface, 'EdgeColor', 'none');
    xlabel('X');
    ylabel('Y');
    zlabel('Surface Height');
    title('Reconstructed Surface from Zernike Coefficients');
    colormap jet; % Apply color map
    colorbar; % Show color bar
    view(3); % Set 3D view
end

function Z = zernike_polynomial(n, m, rho, theta)
    % Function to calculate the Zernike polynomial Z_n^m(rho, theta)
    R = zernike_radial(n, abs(m), rho);
    if m >= 0
        Z = R .* cos(m * theta);
    else
        Z = R .* sin(abs(m) * theta);
    end
end

function R = zernike_radial(n, m, rho)
    % Function to calculate the radial part of the Zernike polynomial
    R = zeros(size(rho));
    for s = 0:((n - m) / 2)
        c = (-1)^s * factorial(n - s) / (factorial(s) * factorial((n + m) / 2 - s) * factorial((n - m) / 2 - s));
        R = R + c * rho.^(n - 2 * s);
    end
end