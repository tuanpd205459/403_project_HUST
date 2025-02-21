%% Xóa dữ liệu trong workspace và đóng tất cả các figure
clear;
close all;
clc;

%% Additional file
filePath = 'C:\path\to\your\folder';
addpath(filePath);

%% Biến toàn cục (tham số đầu vào các hàm)
DPD = 25;
he_so = 1;
poly_order = 3;
%chọn vẽ HCN hay tròn #1: HCN   #0: tròn
chon_vung_HCN = 1;

maxIntensity = 100000;  % Cường độ tối đa để hiển thị ảnh


%% Đọc ảnh hologram đầu vào:
%   inputManual : logical
%       1 - Đọc file bằng tay từ `filePath`
%       0 - Tự động tìm file mới nhất trong `folder_path`
inputManual = 0;
folder_path = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Thi nghiem\6-12';
filePath = '8.bmp';
hologram = processing.loadHologram(inputManual, filePath, folder_path);

wrappedPhase = processing.processFourier(hologram);

% Cách sử dụng   unwrapping.unwrapPhase:
%   unwrapped_Phase = unwrapPhase(wrappedPhase);
%   unwrapped_Phase = unwrapPhase(wrappedPhase, methodGroup);
%   unwrapped_Phase = unwrapPhase(wrappedPhase, methodGroup, methodType);
%
% Đầu vào:
%   - wrappedPhase  : Ma trận pha bị quấn (wrapped phase).
%   - methodGroup   : Nhóm thuật toán tháo gỡ pha (string):
%       + 'poisson'   : Phương pháp Poisson (mặc định dùng DCT).
%       + 'ls'        : Least-Squares (LS).
%       + 'tie'       : Transport of Intensity Equation (TIE).
%       + 'linh'      : Phương pháp của a Linh.
%       + '2dweight'  : Phương pháp 2D weighted phase unwrapping.
%   - methodType    : Phương pháp con (chỉ áp dụng cho LS và TIE) (string):
%       + 'dct'   : Dùng Discrete Cosine Transform (DCT) (mặc định).
%       + 'fft'   : Dùng Fast Fourier Transform (FFT).
%       + 'iter'  : Sử dụng phương pháp có lặp (iterative).
%
% Đầu ra:
%   - unwrapped_Phase: Ma trận pha đã tháo gỡ.
% Ví dụ:
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'ls', 'dct'); % LS với DCT
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'tie', 'fft'); % TIE với FFT
%   unwrapped_Phase = unwrapPhase(wrappedPhase, 'linh'); % Phương pháp của a Linh
%   unwrapped_Phase = unwrapPhase(wrappedPhase, '2dweight'); % 2D weighted phase unwrapping
%
% Nếu không nhập methodGroup và methodType, mặc định sử dụng 'ls' với 'dct'.

methodGroup = 'poisson';
methodType ='';
unwrapped_Phase = unwrapping.unwrapPhase(wrappedPhase, methodGroup);


wavelength = 633e-9;
%he_so = DPD;

% DPD=1; % Do phong dai cua he quang
r = (unwrapped_Phase .* wavelength * he_so) / (4*pi);
r = r * 10^6;  % Initially convert to micrometers

offSet = 10;      
r = r(offSet:end-offSet,offSet:end-offSet);  % cắt/chọn vùng để ve do thi

% r = my_offSet(r);           % xoay ngang mặt phẳng tái tạo
temp_r = r;


imagesc(r);
hold on;                    
positionLine = myDrawLine();  
title('Mặt phẳng pha theo chiều'); % Đặt tiêu đề cho hình ảnh
%colormap(gray);              

%fprintf('\n %.2f, %.2f, %.2f, %.2f\n', positionLine(1,1), positionLine(1,2), positionLine(2,1), positionLine(2,2));
crossLine = myCrossSection(r, positionLine(1,1), positionLine(1,2), positionLine(2,1), positionLine(2,2));
meanLine = myMeanLine(crossLine, poly_order);
[Ra, Ra_line] = myCalcRa(crossLine, meanLine, DPD);
%

numPixels = length(crossLine);
x_real2D = (1:numPixels)*3.45/DPD;


%%
Rz = myCalcRz(crossLine,meanLine);
% fprintf('Độ nhám bề mặt Rz so với đường trung bình: %.4f\n', Rz);
%% Tính Sz
Sz=myCalcSz(r,poly_order);

%% Độ nhám trung bình Sa
Sa = myCalcSa(r, poly_order);
% fprintf('Sa: %.4f %s\n', Sa, dimensional);

%% Độ nhám Sq
Sq = myCalcSq(r, poly_order);
% fprintf('Sq: %.4f %s\n', Sq, dimensional);

%% Chuyển đổi đơn vị 
if max(r(:)) < 1
    r = r * 10^3;  % Convert to nanometers
    dimensional = 'nanomet';
    Ra = Ra * 10^3;
    Rz = Rz * 10^3;
    Sa = Sa * 10^3;
    Sq = Sq * 10^3;
    Sz=Sz*10^3;
    meanLine = meanLine * 10^3;
    Ra_line = Ra_line * 10^3;
    crossLine = crossLine * 10^3;
else
    dimensional = 'micromet';
end

%% Vẽ đồ thị 2D
figure();
plot(x_real2D, crossLine, 'k');
hold on;
plot(x_real2D, meanLine, 'b-', 'LineWidth', 1);
plot(x_real2D, Ra_line, 'r-', 'LineWidth', 1);
legend('Surface', 'Mean Curve', ['Ra: ', num2str(Ra), ' ', dimensional]);
title('Surface Roughness and Mean Curve');
xlabel('Micromet');
ylabel(['Height (', dimensional, ')']);

% Chèn giá trị Ra và Rz vào biểu đồ
Ra_text = ['Ra: ', num2str(Ra, '%.4f'), ' ', dimensional];
Rz_text = ['Rz: ', num2str(Rz, '%.4f'), ' ', dimensional];

% Lấy giới hạn của trục x và trục y
x_limits = xlim;
y_limits = ylim;

% Vị trí hiển thị ở góc dưới bên trái
x_pos = x_limits(1) + 0.05 * (x_limits(2) - x_limits(1));
y_pos_Rz = y_limits(1) + 0.05 * (y_limits(2) - y_limits(1));
y_pos_Ra = y_pos_Rz + 0.05 * (y_limits(2) - y_limits(1)); % Điều chỉnh y_pos_Rz nếu cần thiết

% Chèn văn bản vào biểu đồ
text(x_pos, y_pos_Ra, Ra_text, 'Color', 'k', 'FontSize', 12);
text(x_pos, y_pos_Rz, Rz_text, 'Color', 'k', 'FontSize', 12);

hold off;

%% Tạo figure 3D với đơn vị thích hợp

[N, M] = size(r);
x_min = min(min(r()));
x_max = max(max(r()));
%z_max = max(max(r()));

figure();
x_real = (1:M) * 3.45/DPD; % Tính toán giá trị thực tế của trục x (đơn vị micromet)
y_real = (1:N) * 3.45/DPD; % Tính toán giá trị thực tế của trục y
%x_real = (1:M) ; % Tính toán giá trị thực tế của trục x (đơn vị micromet)
%y_real = (1:N) ; % Tính toán giá trị thực tế của trục y
[X, Y] = meshgrid(x_real, y_real);
% Vẽ biểu đồ mesh
mesh(X, Y, r); 
title('3D Surface'); % Tiêu đề của biểu đồ
colormap(jet);    % Áp dụng bảng màu "jet"
colorbar();
clim([x_min x_max]); %Đặt giới hạn cho cột colorbar
xlabel('x (\mum)'); % Nhãn trục x với đơn vị micromet
ylabel('y (\mum)'); % Nhãn trục y với đơn vị micromet
zlabel(['Height (', dimensional, ')']); % Nhãn trục z với đơn vị thích hợp
hold on;
% my3Dplane(huongMCN, positionLine, DPD);
% hold off;
%%

%%
% Chèn giá trị Sa và Sq vào biểu đồ
Sa_text = ['Sa: ', num2str(Sa, '%.4f'), ' ', dimensional];
Sq_text = ['Sq: ', num2str(Sq, '%.4f'), ' ', dimensional];
Sz_text=['Sz: ', num2str(Sz, '%.4f'), ' ', dimensional];
dim1 = [0.15, 0.9, 0, 0]; % Vị trí và kích thước hộp chứa Sa
dim2 = [0.15, 0.85, 0, 0]; % Vị trí và kích thước hộp chứa Sq
dim3=[0.15, 0.8, 0, 0]; % Vị trí và kích thước hộp chứa Sz
annotation('textbox', dim1, 'String', Sa_text, 'FitBoxToText', 'on', ...
    'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');

annotation('textbox', dim2, 'String', Sq_text, 'FitBoxToText', 'on', ...
    'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');
% annotation('textbox', dim3, 'String', Sz_text, 'FitBoxToText', 'on', ...
%     'EdgeColor', 'none', 'FontSize', 12, 'Color', 'k');

% %%
% z_map = unwrapped_Phase;
% [hang, cot] = size(z_map);
% if (hang > cot)
%     new_z_map = z_map(hang/2-cot/2: hang/2+cot/2 , :);
% else 
%     new_z_map = z_map(: , cot/2-hang/2: cot/2+hang/2);
% end
% figure;
% surf(new_z_map);
% title("new z map");

%% m, n indices
% coeff = zeros(1, 2);
% coeff(1) = 10; coeff(2) = 5;
% [output_coeff, z_recon_map] = ZernikeLegendreFit(new_z_map, "2indices", coeff);
% 
% % figure;
% % subplot(131); imagesc(z_map); colorbar();
% % subplot(132); imagesc(z_recon_map); colorbar();
% % subplot(133); imagesc(z_map - z_recon_map); colorbar();
% 
% % figure;
% % subplot(121); imagesc(output_coeff{1}'); colorbar();
% % subplot(122); imagesc(output_coeff{2}'); colorbar();
% % figure;
% % surf(z_map);
% % title('z map');
% figure;
% surf(z_recon_map);
% title('Tái tạo bề mặt sóng từ 3D');
% xlabel('X');
% ylabel('Y');
% zlabel('Độ lệch pha');
% colormap(jet);    % Áp dụng bảng màu "jet"
% colorbar();

%% Fringe index
% coeff = 100;
% [output_coeff, z_recon_map] = ZernikeLegendreFit(z_map, "fringe", coeff);   
% 
% % figure;
% % subplot(131); imagesc(z_map); colorbar()
% % subplot(132); imagesc(z_recon_map); colorbar();
% % subplot(133); imagesc(z_map - z_recon_map); colorbar();
% % figure;
% % plot(output_coeff{1});
% 
% fprintf('He so: %.1f \n', output_coeff{1});
% 
% 
% %% Numerical reconstruction (In fact, no need because z0=0 as we record the hologram at the focus plane)
% h = 0.633e-6;  % Bước sóng ánh sáng
% k = 2*pi/h;    % Số sóng
% z0 = 10;       % Khoảng cách truyền ánh sáng
% 
% [numRows, numCols] = size(r);  % Kích thước của ma trận hình ảnh 'r'
% 
% % Tính toán tọa độ không gian của ảnh (cần xác định các biến 'lengthOx' và 'lengthOy')
% n_pixels = 1:numCols;  % Chỉ số pixel theo chiều ngang
% m_pixels = 1:numRows;  % Chỉ số pixel theo chiều dọc
% 
% % Cần chắc chắn rằng 'lengthOx' và 'lengthOy' là các giá trị chiều dài vật lý của vùng ảnh
% x_coords = h * (-numCols/(2*lengthOx) + (n_pixels-1)/lengthOx);  % Tọa độ x
% y_coords = h * (-numRows/(2*lengthOy) + (m_pixels-1)/lengthOy);  % Tọa độ y
% 
% [xx, yy] = meshgrid(x_coords, y_coords);  % Lưới tọa độ không gian
% 
% % Kiểm tra điều kiện sqrt(1-xx.^2-yy.^2) để tránh giá trị phức
% trans = exp(1i*k*z0*sqrt(1 - xx.^2 - yy.^2));  % Phương pháp tần số góc
% result = r .* trans;  % Nhân với ma trận truyền
% 
% Uf = ifft2(fftshift(result));  % Biến đổi ngược Fourier
% If = Uf .* conj(Uf);  % Cường độ

% % Hiển thị kết quả tái tạo
% figure;
% surf(angle(Uf));  % Hiển thị pha của ảnh
% title("Ảnh 3D sau truyền ngược");
% shading interp;  % Làm mịn bề mặt
% lighting phong;  % Shading phong để ánh sáng tốt hơn
% camlight headlight;  % Thêm nguồn sáng từ góc nhìn của camera
% 
% figure;
% imshow(If, []);  % Hiển thị cường độ của kết quả tái tạo
% colormap(gray);  % Áp dụng bản đồ màu xám
% title('Final result intensity');


%% Xoay mặt phẳng tái tạo nằm ngang
function offSet_surface = my_offSet(r)
    % Lấy cột giữa của ma trận r
%   middle = r(:, round(end/2));
  middle = r(round(end/2),:);
    % Trừ cột giữa từ mỗi cột trong ma trận r
    % MATLAB sẽ tự động lặp lại middle_column cho mỗi cột của r
    offSet_surface = r - middle;
end





%% Hàm unwrap của PCA 
%{
function up= uphase(wp)
    n=zeros(size(wp));
    for i=2:length(wp)
        n(i)=floor((wp(i)-wp(i-1))/(2*pi)+0.5)+n(i-1);   
    end
     up=-2*pi*n+wp;
end
%}

%% function tính Mean line
function outMeanLine = myMeanLine(inputRow, poly_order)
    %global poly_order
    N = length(inputRow);
    x = 1:N;
    p = polyfit(x, inputRow, poly_order);
    % Tính toán giá trị của đường trung bình cong
    outMeanLine = polyval(p, x);
end

%% Hàm chuyển sang grayscale
function output = myConvGrayScale(inputImage)
    if size(inputImage, 3) > 1
        inputImage = rgb2gray(inputImage);
    end
    output = double(inputImage); 
end
%% Hàm tính Ra
function [Ra, Ra_line] = myCalcRa(inputLine, input_meanline, DPD)
    if isrow(inputLine)
        inputLine = inputLine'; 
    end
    if isrow(input_meanline)
        input_meanline = input_meanline';  
    end
    lr = length(inputLine) * 3.45/DPD;  % Length of the surface in micrometers
    u = linspace(1, lr, length(inputLine));  
    
    % tính Ra = tích phân hình thang
    Ra = (1/lr) * trapz(u, abs(inputLine - input_meanline));
    
    % Kết quả
  %  disp(['Surface Roughness Ra: ', num2str(Ra), ' micrometers']);
    
    Ra_line = Ra + input_meanline;
end
%% Hàm tính Rz
function output = myCalcRz(inputSuface, inputMeanLine)
    differenceValue = inputSuface - inputMeanLine;
    sorted_differenceValue = sort(differenceValue(:));
    
    top5_sum = sum(sorted_differenceValue(end-4:end));
    bottom5_sum = sum(sorted_differenceValue(1:5));
    output = (top5_sum - bottom5_sum) / 5;

end

%% Hàm vẽ đường thẳng
function positionLine = myDrawLine()
    roiLine = drawline('Color','r');
    wait(roiLine);
    positionLine = round(roiLine.Position);
end
%% hàm vẽ hình chữ nhật
function [pos, xRec, yRec, widthRec, heightRec] = myDrawRec()
%     roi = drawrectangle();  % Rectangle

% Vẽ hình chữ nhật để chọn ROI
roi = drawrectangle();

% Lấy tọa độ tâm ban đầu
centerRec = round([roi.Position(1) + roi.Position(3)/2, roi.Position(2) + roi.Position(4)/2]);

% Vẽ dấu cộng tại tâm
hold on;
centerMarker = plot(centerRec(1), centerRec(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% Thiết lập callback để cập nhật tâm trong quá trình di chuyển hình chữ nhật
addlistener(roi, 'MovingROI', @(src, evt) updateCenterRectangle(src, centerMarker));




    wait(roi);  % Double-click to confirm the ROI
    
    % Get the position of the rectangle [x, y, width, height]
    pos = round(roi.Position);
    xRec = pos(1);
    yRec = pos(2);
    widthRec = pos(3);
    heightRec = pos(4);
end


% Hàm cập nhật tâm khi di chuyển ROI
function updateCenterRectangle(roi, centerMarker)
    % Cập nhật vị trí của dấu cộng theo tọa độ tâm mới
    centerX = roi.Position(1) + roi.Position(3)/2;
    centerY = roi.Position(2) + roi.Position(4)/2;
    centerMarker.XData = centerX;
    centerMarker.YData = centerY;
    drawnow;
end

%% hàm vẽ đường tròn
function [centerCir_X, centerCir_Y, radiusCir] = myDrawCir()
    %roi = drawcircle();

    % Vẽ vòng tròn để chọn ROI
    roi = drawcircle();
    
    % Lấy tọa độ ban đầu của tâm và bán kính
    centerCir = round(roi.Center);
    
    % Vẽ dấu cộng tại tâm
    hold on;
    centerMarker = plot(centerCir(1), centerCir(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    
    hold off;
    
    % Thiết lập callback để cập nhật tâm trong lúc di chuyển vòng tròn
    addlistener(roi, 'MovingROI', @(src, evt) updateCenter(src, centerMarker));
    
    % Hàm cập nhật tâm trong quá trình di chuyển ROI



    wait(roi);
    centerCir_X = round(roi.Center(1));  % Tọa độ x của trung tâm hình tròn
    centerCir_Y = round(roi.Center(2));  % Tọa độ y của trung tâm hình tròn
    radiusCir = round(roi.Radius);   
end

    function updateCenter(roi, centerMarker)
        % Cập nhật vị trí của dấu cộng theo tọa độ tâm mới
        centerMarker.XData = roi.Center(1);
        centerMarker.YData = roi.Center(2);
        drawnow;
    end



%% Hàm mặt cắt ngang
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
%% Hàm tính sai lệch z(x,y) - để tính các thông số nhám 3D
function difference_Matrix = myDifference (inputMatrix, poly_order)
    % Khởi tạo mảng để lưu kết quả Sz cho từng hàng
    difference_Matrix = zeros(size(inputMatrix));
    for i = 1:size(inputMatrix, 1)
        row_vector = inputMatrix(i, :);
        % Tính giá trị trung binh của hàng i
        meanLine_row = myMeanLine(row_vector,poly_order);
        % Tính sai lech hang i
        difference_Matrix (i,:) = inputMatrix(i,:) - meanLine_row;
    end
end

%% Hàm tính Sa
function Sa = myCalcSa (inputMatrix, poly_order)
    [numRows, numCols] = size(inputMatrix);
    z_xy = myDifference(inputMatrix, poly_order);   % Tính sai lệch các giá trị
    abs_z = abs(z_xy(:)); 
    Sa = sum(abs_z(:)) / (numRows * numCols); 
    
end

%% Hàm tính Sq
function Sq = myCalcSq (inputMatrix, poly_order)
    [numRows, numCols] = size(inputMatrix);
    z_xy = myDifference(inputMatrix, poly_order);   % Tính sai lệch các giá trị
    abs_z = abs(z_xy(:)); 
    Sq = sqrt(1/(numRows*numCols) * sum(abs_z(:).^2));
    
end

%% Hàm tính Sz
function Sz_avg = myCalcSz (inputMatrix, poly_order)
    % Khởi tạo mảng để lưu kết quả Sz cho từng hàng
    Sz_row = zeros(size(inputMatrix, 1), 1);
    
    for i = 1:size(inputMatrix, 1)
        row_vector = inputMatrix(i, :);
        
        % Tính giá trị trung bình của hàng i
        meanLine_row = myMeanLine(row_vector, poly_order);
        
        % Tính độ lệch bình phương
        Rz_row = myCalcRz(row_vector, meanLine_row);
        
        % Tính tổng và căn bậc hai của tổng cho hàng i
        sum_Rz_row = sum(Rz_row);
        Sz_row(i) = sqrt(sum_Rz_row);
    end
    
    % Tính giá trị trung bình của Sz cho các hàng
    Sz_avg = mean(Sz_row);
    
    disp(['Sz theo các hàng: ', num2str(Sz_avg)]);
end
%% Hàm vẽ mặt cắt ngang 3D
function my3Dplane(input3D, inputPosition, DPD)
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


 %% Bắt đầu thuật toán PCA
function wrappedPhase = performPCA(roiContent, fourierTransform)
    % Lấy mẫu và tính toán phần mũ
    compensationSpectrum = roiContent;
    [heightRec,widthRec] = size(roiContent);
    % IFFT the cropped spectrum to get a subsample hologram
    subsampledHologram = ifft2(fftshift(compensationSpectrum));
    
    rawSubHolo = subsampledHologram;
    
    % Get the exponential term
    expTerm = exp(-1i * angle(subsampledHologram));
    
    %% SVD 
    % SVD and get the first dominant principal component
    [U,S,V] = svd(expTerm);
    %Y = U*S*V';
    
    [m,n]=size(S); 
    SS=zeros(m,n); 
    SS(1,1)=S(1,1);
    Z= U*SS*V';
    SS(2,2)=S(2,2);
    Z2= U*SS*V';
    SS(3,3)=S(3,3);
    SS(4,4)=S(4,4);
    SS(5,5)=S(5,5);
    SS(6,6)=S(6,6);
%    SS(7,6)=S(7,7);
    SS(7,7)=S(7,7);
    Z3= U*SS*V';
    
%     figure;
%     subplot(1,3,1),imshow(angle(Z),[]); 
%     colormap(gray);
%     title('1st PC');
%     subplot(1,3,2),imshow(angle(Z2),[]);
%     colormap(gray);
%     title('1st+2nd PCs');
%     subplot(1,3,3),imshow(angle(Z3),[]);
%     colormap(gray);
%     title('1st-7th PCs');
    
    
    % Least-squares fitting to the two dominants singular vector
    UNU1 = unwrap(angle(U(:,1)));
    A = polyfit(1:heightRec,UNU1',2);
    NewUN1=polyval(A,1:heightRec);
    NEWU1=exp(1i.*(NewUN1'));
    U(:,1)=NEWU1;
    
    UNV1 = unwrap(angle(V(:,1)));
    A = polyfit(1:widthRec,UNV1',2);
    NewUNV1 = polyval(A,1:widthRec);
    NEWV1 = exp(1i.*((NewUNV1')));
    V(:,1) = NEWV1;
    
    % The the aberration term
    SS=zeros(m,n); 
    SS(1,1)=S(1,1);
    Z= U*SS*V';
    
%     figure
%     subplot(1,3,1),
%     imshow(angle(rawSubHolo),[]);
%     title('Raw subsampled phase');
%     subplot(1,3,2)
%     imshow(angle(Z),[]);
%     title('Conjugated phase extracted');
  %  outPutPCA = Z;

    %% Bù đắp bằng cách sử dụng thuật ngữ sai lệch trích xuất 
    % (trung tâm phổ được thực hiện sau khi bù đắp PCA, 
    % trong khi trong bài báo được thực hiện trước đó), 
    % thực tế không ảnh hưởng đến kết quả.
    [numRows, numCols] = size(fourierTransform);  % Lấy kích thước ảnh Fourier
    
    % Tính toán vị trí trung tâm mới
    centerX = round((numRows - heightRec) / 2);  % Tọa độ x của tâm
    centerY = round((numCols - widthRec) / 2); % Tọa độ y của tâm
    
    newFourierTransform = zeros(size(fourierTransform));
    
    % Multiple the conjugation of aberration term with the subsample hologram
    % FFT and replace the corresponding orginal region of the spectrum
    newFourierTransform(centerX:centerX + heightRec - 1, centerY:centerY + widthRec - 1) = ...
        fftshift(fft2(((rawSubHolo) .* Z)));
    
    
    subplot(1, 3, 3);
    imshow(log(abs(newFourierTransform)+1), []);
    title('Compensated Spectrum');
    
    % IFFT to retrieve phase
    finalPhase = ifft2(ifftshift(newFourierTransform));  
    wrappedPhase = angle(finalPhase);  
    figure;
    mesh(wrappedPhase); title('Wrapped Phase');

end




function finds_center_square(big_square_size, small_square_size,step)
%   Kích thước hình vuông lớn và nhỏ
%   big_square_size = 10;  % Hình vuông lớn có kích thước 10x10
%   small_square_size = 6; % Hình vuông nhỏ có kích thước 4x4

% Vẽ hình vuông lớn
figure; % Tạo một cửa sổ mới
hold on;

% Vẽ hình vuông lớn
rectangle('Position', [0 0 big_square_size big_square_size], 'EdgeColor', 'b', 'LineWidth', 2);

% Duyệt qua tất cả các vị trí có thể để đặt hình vuông nhỏ 4x4 trong hình vuông lớn
for x_start = 0:step:(big_square_size - small_square_size)
    for y_start = 0:step:(big_square_size - small_square_size)
        % Vẽ hình vuông nhỏ tại mỗi vị trí
     %  rectangle('Position', [x_start y_start small_square_size small_square_size], 'EdgeColor', 'r', 'LineWidth', 1);
        
        % Tính tọa độ tâm của hình vuông nhỏ
        x_center = x_start + small_square_size/2;
        y_center = y_start + small_square_size/2;
        
        % Hiển thị tâm bằng dấu chấm tròn màu đỏ
        plot(x_center, y_center, 'ro', 'MarkerFaceColor', 'r');
    end
end

% Đặt giới hạn cho trục
axis([0 big_square_size 0 big_square_size]);
axis square;
grid on;

title('Tất cả các hình vuông 4x4 và tâm của chúng trong hình vuông lớn 10x10');
xlabel('Trục X');
ylabel('Trục Y');
hold off;
end

function surface = reconstruct_wavefront(zernike_coeffs, order, grid_size)
    % Hàm tái tạo bề mặt sóng từ các hệ số đa thức Zernike
    % zernike_coeffs: Hệ số đa thức Zernike
    % order: Bậc của các đa thức Zernike
    % grid_size: Kích thước của lưới điểm trên mặt phẳng 

    % Khởi tạo mặt phẳng lưới
    [X, Y] = meshgrid(linspace(-1, 1, grid_size), linspace(-1, 1, grid_size));
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X);

    % Khởi tạo bề mặt sóng
    surface = zeros(size(X));

    % Lặp qua các bậc và hệ số để tái tạo bề mặt sóng
    index = 1;
    for n = 0:order
        for m = -n:2:n
            if index <= length(zernike_coeffs)
                % Lấy hệ số tương ứng
                Z = zernike_coeffs(index);
                
                % Tính giá trị Zernike polynomial
                Zmn = ZernikePolynomial(n, m, R, Theta);
                
                % Cộng dồn vào bề mặt sóng
              %  surface = surface + Z * Zmn;
                surface = surface + Z * Zmn;

                index = index + 1;
            end
        end
    end
end

function Zmn = ZernikePolynomial(n, m, R, Theta)
    % Hàm tính Zernike polynomial
    % n: Bậc n của Zernike
    % m: Bậc m của Zernike
    % R: Bán kính
    % Theta: Góc theta
    
    % Tính radial Zernike polynomial
    RadialPoly = zeros(size(R));
    for k = 0:floor((n-abs(m))/2)
        c = (-1)^k * factorial(n-k) / (factorial(k) * factorial((n + abs(m))/2 - k) * factorial((n - abs(m))/2 - k));
        RadialPoly = RadialPoly + c * R.^(n - 2*k);
    end
    
    % Tính Zernike polynomial
    if m >= 0
        Zmn = RadialPoly .* cos(m * Theta);
    else
        Zmn = RadialPoly .* sin(abs(m) * Theta);
    end
end

%% Tính toán đa thức Zernike
function radial = zernike_radial(r,n,m)
    % Functions required for use: elliptical_crop
%     hàm tính toán đa thức Zernike
%     Đầu vào
%         r:      bán kính
%         n:      bậc quang sai
%         m:      bậc phương vị
%     Đầu ra:
%         giá trị Đa thức Zernike

    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0 || m < 0
        error('n and m must both be positive in radial function')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end
    if n == m
        radial = r.^n;
    elseif n - m == 2
        radial = n*zernike_radial(r,n,n)-(n-1)*zernike_radial(r,n-2,n-2);
    else
        H3 = (-4*((m+4)-2)*((m+4)-3)) / ((n+(m+4)-2)*(n-(m+4)+4));
        H2 = (H3*(n+(m+4))*(n-(m+4)+2)) / (4*((m+4)-1))  +  ((m+4)-2);
        H1 = ((m+4)*((m+4)-1) / 2)  -  (m+4)*H2  +  (H3*(n+(m+4)+2)*(n-(m+4))) / (8);
        radial = H1*zernike_radial(r,n,m+4) + (H2+H3 ./ r.^2).*zernike_radial(r,n,m+2);
        
        % Fill in NaN values that may have resulted from DIV/0 in prior
        % line. Evaluate these points directly (non-recursively) as they
        % are scarce if present.
        
        if sum(sum(isnan(radial))) > 0
            [row, col] = find(isnan(radial));
            c=1;
            while c<=length(row)
                x = 0;
                for k = 0:(n-m)/2
                    ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                    x = x + ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                end
                radial(row(c),col(c)) = x;
                c=c+1;
            end
        end

    end

end
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


