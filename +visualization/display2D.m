function display2D(surfaceParams)

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
    
    % Số pixel trên mặt cắt
    numPixels = length(crossLine);
    x_real2D = (1:numPixels) * 3.45 / DPD; % Trục x theo micromet

    % Vẽ đồ thị 2D
    figure();
    plot(x_real2D, crossLine, 'k'); hold on;
    plot(x_real2D, meanLine, 'b-', 'LineWidth', 1);
    plot(x_real2D, Ra_line, 'r-', 'LineWidth', 1);
    
    % Thêm chú thích
    legend('Surface', 'Mean Curve', ['Ra: ', num2str(Ra), ' ', dimensional], 'Location', 'best');
    title('Độ nhám bề mặt và đường trung bình');
    xlabel('Micromet');
    ylabel(['Height (', dimensional, ')']);

    % Lấy giới hạn của trục x và y
    x_limits = xlim;
    y_limits = ylim;

    % Tính vị trí hiển thị văn bản
    x_pos = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1));
    y_pos_Rz = y_limits(1) + 0.05 * (y_limits(2) - y_limits(1));
    y_pos_Ra = y_pos_Rz + 0.05 * (y_limits(2) - y_limits(1)); % Điều chỉnh khoảng cách

    % Chèn giá trị Ra và Rz vào biểu đồ
    text(x_pos, y_pos_Ra, sprintf('Ra: %.4f %s', Ra, dimensional), 'Color', 'k', 'FontSize', 12);
    text(x_pos, y_pos_Rz, sprintf('Rz: %.4f %s', Rz, dimensional), 'Color', 'k', 'FontSize', 12);

    hold off;
end
