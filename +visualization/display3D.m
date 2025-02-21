function display3D(surfaceParams, enableDisplayCrossSection3D)
    % Kiểm tra đầu vào hợp lệ
    if isempty(enableDisplayCrossSection3D)
        enableDisplayCrossSection3D = false;
    end
%     surfaceParams = {reconSurface, Ra, Rz, Sa, Sq, Sz, Ra_line, positionLine,...
%                                     crossLine, meanLine, dimensional, DPD}; 
    reconSurface = surfaceParams{1};
    Ra = surfaceParams{2};
    Rz = surfaceParams{3};
    Sa = surfaceParams{4};
    Sq = surfaceParams{5};
    Sz = surfaceParams{6};
    Ra_line = surfaceParams{7};
    positionLine = surfaceParams{8};
    crossLine = surfaceParams{9};
    meanLine = surfaceParams{10};
    dimensional = surfaceParams{11};
    DPD = surfaceParams{12};

    % Kích thước của bề mặt tái tạo
    [N, M] = size(reconSurface);
    
    % Lấy giá trị min/max cho trục màu
    z_min = min(reconSurface(:));
    z_max = max(reconSurface(:));

    % Chuyển đổi tọa độ thực tế
    x_real = (1:M) * 3.45 / DPD; % Trục x (micromet)
    y_real = (1:N) * 3.45 / DPD; % Trục y (micromet)
    [X, Y] = meshgrid(x_real, y_real);

    % Vẽ đồ thị 3D
    figure();
    mesh(X, Y, reconSurface);
    title('3D Surface'); 
    colormap(jet);
    colorbar();
    clim([z_min, z_max]); % Đặt giới hạn cho colorbar
    xlabel('x (\mum)'); 
    ylabel('y (\mum)');
    zlabel(['Height (', dimensional, ')']);

    hold on;

    %% Chèn giá trị Sa, Sq, Sz vào biểu đồ
    Sa_text = sprintf('Sa: %.4f %s', Sa, dimensional);
    Sq_text = sprintf('Sq: %.4f %s', Sq, dimensional);
    Sz_text = sprintf('Sz: %.4f %s', Sz, dimensional);

    % Vị trí annotation
    dim1 = [0.15, 0.88, 0, 0]; % Sa
    dim2 = [0.15, 0.83, 0, 0]; % Sq
    dim3 = [0.15, 0.78, 0, 0]; % Sz

    % Thêm annotation vào hình
    annotation('textbox', dim1, 'String', Sa_text, 'FitBoxToText', 'on', ...
        'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');
    annotation('textbox', dim2, 'String', Sq_text, 'FitBoxToText', 'on', ...
        'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');
    annotation('textbox', dim3, 'String', Sz_text, 'FitBoxToText', 'on', ...
        'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');

    %% chèn mặt cắt ngang vào biểu đồ
    if enableDisplayCrossSection3D
        visualization.displayCrossSection3D(reconSurface, positionLine, DPD);
    end
    hold off;
end
