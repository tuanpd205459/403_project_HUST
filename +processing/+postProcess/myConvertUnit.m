function [reconSurface, dimensional] = myConvertUnit(reconSurface)
    % Kiểm tra giá trị lớn nhất để xác định đơn vị
%     if (max(reconSurface(:))-min(reconSurface(:))) < 10
    if max(abs(reconSurface(:))) < 1
        scaleFactor = 10^3; % Chuyển từ micromet sang nanomet
        dimensional = 'nanomet';
    else
        scaleFactor = 1;  
        dimensional = 'micromet';
    end

    % Chuyển đổi đơn vị
    reconSurface = reconSurface * scaleFactor;
end
