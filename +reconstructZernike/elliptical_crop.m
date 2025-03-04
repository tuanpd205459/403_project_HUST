function cropped_im = elliptical_crop(im,crop_frac)
    % Hàm crop ảnh tương ứng đường tròn
    
    if crop_frac < 0 || crop_frac > 1
        error('crop_frac must have value between 0 and 1')
    end

    cropped_im = im;
    center_x = (size(im,2)+1)/2;
    center_y = (size(im,1)+1)/2;
    radius_x = (size(im,2)-center_x)*crop_frac;
    radius_y = (size(im,1)-center_y)*crop_frac;

    for row = 1:size(im,1)
        for col = 1:size(im,2)
            if sqrt((row-center_y)^2/radius_y^2 + (col-center_x)^2/radius_x^2) > 1
                cropped_im(row,col) = nan;
            end
        end
    end
    
    % Necessary because of potential DIV/0 behavior
    if radius_x == 0 || radius_y == 0
        cropped_im = nan(size(cropped_im));
    end

end