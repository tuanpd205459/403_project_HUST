%% Hàm vẽ mặt cắt ngang 3D
function displayCrossSection3D(input3D, inputPosition, DPD)
    z_max = max(max(input3D()));
    z_min = min(min(input3D()));

    % Tính toán giá trị thực tế của trục x (đơn vị micromet)
    x1 = inputPosition(1,1) * 3.45/DPD;
    y1 = inputPosition(1,2) * 3.45/DPD;
    x2 = inputPosition(2,1) * 3.45/DPD;
    y2 = inputPosition(2,2) * 3.45/DPD;
    
    normal = cross([x2- x1, y2-y1, 0], [0, 0, z_max]);
    % phương trình ax1 + by1 + cz1 = d , di qua diem A(x1, y1, 0)
    a = normal(1);
    b = normal(2);
    %c = normal(3);
    d = a*x1 + b*y1;

    % Tạo lưới điểm cho mặt phẳng
    
    [xp, zp] = meshgrid(linspace(x1, x2, 100), linspace(z_min, z_max, 100));
    % Tính toán các giá trị y từ phương trình mặt phẳng ax1 + by1 = d
    yp = (d - a*xp) / b; 
    
    % Vẽ mặt phẳng song song với trục Oz
    surf(xp, yp, zp, 'FaceColor', 'blue', 'EdgeColor', 'black');

    % Tắt lưới cho mặt phẳng
    %shading interp; % hoặc shading flat

end