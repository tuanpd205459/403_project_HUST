clc; clear; close all

x = sin(linspace(0, 6*pi, 201));
y = sin(linspace(0, 3.5*pi, 198))/2 + 0.3;
z_map = y' * x;

%% m, n indices
coeff = zeros(1, 2);
coeff(1) = 20; coeff(2) = 10;
[output_coeff, z_recon_map] = ZernikeLegendreFit(z_map, "2indices", coeff);

figure;
subplot(131); surf(z_map); colorbar()
subplot(132); surf(z_recon_map); colorbar();
subplot(133); surf(z_map - z_recon_map); colorbar();
figure;
subplot(121); surf(output_coeff{1}'); colorbar();
disp(output_coeff{1}');


