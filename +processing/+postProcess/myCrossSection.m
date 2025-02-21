function crossLine = myCrossSection(inputSurface, x1, y1, x2, y2)
    if(x1 ~= x2)
        num_samples = abs(x2- x1);
        % Tạo các giá trị x và y dọc theo đường thẳng sử dụng nội suy
        x = linspace(x1, x2, num_samples);
        y = linspace(y1, y2,num_samples);
        
        % Nội suy các giá trị cường độ dọc theo đường thẳng
        crossLine = interp2(inputSurface, x, y);
%         Vẽ mặt cắt ngang của cường độ
        
        x_micromet = (1:num_samples)*3.45;
        figure;
        plot(x_micromet, crossLine);
        title('MCN pha');
        xlabel('x \mum');
        ylabel('y (nanomet)');
        
    else
        num_samples = abs(y2- y1);
        % Tạo các giá trị x và y dọc theo đường thẳng sử dụng nội suy
        x = linspace(x1, x2, num_samples);
        y = linspace(y1, y2,num_samples);
        
        % Nội suy các giá trị cường độ dọc theo đường thẳng
    %    crossLine = interp2(inputSurface, y, x); 
         crossLine = interp2(inputSurface, x, y);
    %     crossLine = crossLine([2,1]);
    end
end