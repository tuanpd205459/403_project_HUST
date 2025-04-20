%% start

%% Xóa dữ liệu trong workspace và đóng tất cả các figure
clear;
close all;
clc;

%% Additional file (thêm file ảnh cần thiết - nếu nằm ngoài thư mục source code)
filePath = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Ảnh mẫu';
addpath(filePath);

%% Thêm đường dẫn tới file chứa hệ số Zernike


%% Biến toàn cục (tham số đầu vào các hàm)
DPD = 25;
he_so = 1;
poly_order = 3;
maxIntensity = 100000 ;  % Cường độ tối đa để hiển thị ảnh

%% Đọc ảnh hologram đầu vào:
%
% loadImgHologram - Đọc ảnh hologram từ file hoặc thư mục
%
% Syntax:
%   loadImgHologram(inputManual, filePath, folder_path)
%
% Description:
%   Hàm này đọc ảnh hologram từ đường dẫn cụ thể hoặc tự động lấy ảnh mới nhất trong thư mục.
%
% Input Arguments:
%   inputManual : logical
%       1 - Đọc file bằng tay từ `filePath`
%       0 - Tự động tìm file mới nhất trong `folder_path`
%
%   filePath : string
%       Đường dẫn đầy đủ đến file ảnh (chỉ dùng nếu input file bằng tay`).
%
%   folder_path : string
%       Đường dẫn đến thư mục chứa ảnh (chỉ dùng nếu input file tự động`).
%
% Example:
%   % Đọc ảnh từ đường dẫn cụ thể
%   loadImgHologram(1, 'C:\data\hologram.bmp', '')
%
%   % Tự động lấy ảnh mới nhất từ thư mục
%   loadImgHologram(0, '', 'C:\data')

%   inputManual : logical
%       1 - Đọc file bằng tay từ `filePath`
%       0 - Tự động tìm file mới nhất trong `folder_path`
inputManual = 1;
folder_path = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Thi nghiem\6-12';
% filePath = '8.bmp';

filePath = '41.bmp';

hologram = processing.loadHologram(inputManual, filePath, folder_path);

%% processFourier: thực hiện các phép biến đổi Fourier và chọn vùng bậc +1
% wrappedPhase = processing.processFourier(hologram);
wrappedPhase = processing.processFourier(hologram);



%% Chọn thuật toán unwrapping
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



%%
%he_so = DPD;
wavelength = 633;
% DPD=1; % Do phong dai cua he quang
reconSurface = (unwrapped_Phase .* wavelength .* he_so) / (4*pi);
% reconSurface = reconSurface * 10^6;     
% reconSurface = unwrapped_Phase;
offSet = 10;      
reconSurface = reconSurface(offSet:end-offSet,offSet:end-offSet);  % cắt/chọn vùng để vẽ đồ thị

temp_r = reconSurface;
%%
Z = reconSurface;
% % Kích thước bề mặt
[rows, cols] = size(Z);
% 
% % Tạo lưới tọa độ tương ứng
[X, Y] = meshgrid(1:cols, 1:rows);
% 
% % Vector hóa dữ liệu để phù hợp với fit đa biến
x = X(:);
y = Y(:);
% z = Z(:);
% 
% % Fit mặt phẳng Z = a*x + b*y + c bằng hồi quy tuyến tính
% A = [x, y, ones(size(x))];
% coeff = A \ z;  % giải phương trình bình phương tối thiểu
% 
% % Tính mặt phẳng nghiêng đã khớp
% Zfit = reshape(A * coeff, size(Z));
% 
% % Trừ mặt phẳng nghiêng khỏi bề mặt ban đầu
% Zcorr = Z - Zfit;
Z = imgaussfilt(Z, 1);
% 2. Fit mặt phẳng: z = ax + by + c
A = [X(:), Y(:), ones(numel(X),1)];
coeffs = A \ Z(:);
a = coeffs(1); b = coeffs(2); c = coeffs(3);

% 3. Tính độ lệch vuông góc đến mặt phẳng
denom = sqrt(a^2 + b^2 + 1);
Zcorr = (Z - (a*X + b*Y + c)) / denom;

reconSurface = Zcorr;

% load('averageMatrix.mat');
save("main1.mat");
run('main_tai_tao_pha_bu.m'); 
run('catanh.m');
run('bu_pha.m');
run('main1_zernike.m');
run('bu_pha_zernike.m');
load("bu_pha.mat");
load("bu_pha_zernike.mat");



%% Chuyển đổi đơn vị (sang nanomet nếu cần)
[reconSurface, dimensional] = processing.postProcess.myConvertUnit(reconSurface);

%% Xác định chiều vân (ngang/dọc)
detectFringe = processing.postProcess.detectFringeSobel(reconSurface);
disp(detectFringe);
if strcmpi(detectFringe, 'vân ngang')
    reconSurface = rot90(reconSurface); % Xoay 90 độ theo chiều dương
end
if ~isequal(size(reconSurface), size(Zcorr))
    Zcorr = rot90(Zcorr);
end

%% Post Processing
figure;
imagesc(reconSurface);
hold on;    
title('Mặt phẳng pha'); % Đặt tiêu đề cho hình ảnh

% Vẽ đường thằng cắt ngang
positionLine = processing.postProcess.myDrawLine();  

crossLine = processing.postProcess.myCrossSection(reconSurface, positionLine);

% mcn_da_xoay = processing.postProcess.myCrossSection(Zcorr, positionLine);
%%
% % Vẽ đồ thị 2D
%  numPixels = length(mcn_da_xoay);
%     x_real2D = (1:numPixels) * 3.45 / DPD; % Trục x theo micromet
%     figure();
%     plot(mcn_da_xoay, 'k');
%     title("Mat cat ngang da xoay");

% tinh đường trung bình mean_line
meanLine = processing.postProcess.myMeanLine(crossLine, poly_order);

%% Tính toán độ nhám 2D
[Ra, Ra_line] = roughness.myCalcRa(crossLine, meanLine, DPD);

Rz = roughness.myCalcRz(crossLine, meanLine);

%% Tính toán độ nhám 3D
Sz=  roughness.myCalcSz(reconSurface,poly_order);
% Độ nhám trung bình Sa
Sa = roughness.myCalcSa(reconSurface, poly_order);
% Độ nhám Sq
Sq = roughness.myCalcSq(reconSurface, poly_order);

%% Tổng hợp thông số độ nhám
surfaceParams = {reconSurface, Ra, Rz, Sa, Sq, Sz, Ra_line, positionLine,...
                                    crossLine, meanLine, dimensional, DPD}; 

%% Vẽ đồ thị 2D
visualization.display2D(surfaceParams);

%% Tạo figure 3D với đơn vị thích hợp
enableDisplayCrossSection3D = false;     %bật tắt mặt cắt ngang ảnh 3D

visualization.display3D(surfaceParams, enableDisplayCrossSection3D);


%% end1

%% Post Processing 2
figure;
imagesc(Zcorr);
hold on;    
title('Mặt phẳng pha'); % Đặt tiêu đề cho hình ảnh

% Vẽ đường thằng cắt ngang
positionLine = processing.postProcess.myDrawLine();  

crossLine = processing.postProcess.myCrossSection(Zcorr, positionLine);

mcn_da_xoay = processing.postProcess.myCrossSection(Zcorr, positionLine);
%%
% Vẽ đồ thị 2D
 numPixels = length(mcn_da_xoay);
    x_real2D = (1:numPixels) * 3.45 / DPD; % Trục x theo micromet
    figure();
    plot(mcn_da_xoay, 'k');
    title("Mat cat ngang da xoay");

% tinh đường trung bình mean_line
meanLine = processing.postProcess.myMeanLine(crossLine, poly_order);

%% Tính toán độ nhám 2D
[Ra, Ra_line] = roughness.myCalcRa(crossLine, meanLine, DPD);

Rz = roughness.myCalcRz(crossLine, meanLine);

%% Tính toán độ nhám 3D
Sz=  roughness.myCalcSz(Zcorr,poly_order);
% Độ nhám trung bình Sa
Sa = roughness.myCalcSa(Zcorr, poly_order);
% Độ nhám Sq
Sq = roughness.myCalcSq(Zcorr, poly_order);

%% Tổng hợp thông số độ nhám
surfaceParams = {Zcorr, Ra, Rz, Sa, Sq, Sz, Ra_line, positionLine,...
                                    crossLine, meanLine, dimensional, DPD}; 

%% Vẽ đồ thị 2D
visualization.display2D(surfaceParams);

%% Tạo figure 3D với đơn vị thích hợp
enableDisplayCrossSection3D = false;     %bật tắt mặt cắt ngang ảnh 3D

visualization.display3D(surfaceParams, enableDisplayCrossSection3D);

%%
% Hàm xử lý dữ liệu từ thư mục
function averageMatrix = processCSVData(folderPath, rowRange, columnRange, numRows, numCols)
    csvFiles = dir(fullfile(folderPath, '*.csv'));
    numFiles = length(csvFiles);
    
    if numFiles == 0
        error('Không tìm thấy file CSV nào trong thư mục: %s', folderPath);
    end
    
    allData = cell(1, numFiles);
    for i = 1:numFiles
        data = readmatrix(fullfile(folderPath, csvFiles(i).name));
        
        if size(data, 1) >= rowRange(2) && size(data, 2) >= columnRange(2)
            allData{i} = data(rowRange(1):rowRange(2), columnRange(1):columnRange(2));
        else
            warning('File %s không đủ kích thước.', csvFiles(i).name);
            allData{i} = NaN(numRows, numCols);
        end
    end
    
    % Tính trung bình cộng
    averageMatrix = mean(cat(3, allData{:}), 3, 'omitnan');
end

function [reconSurface, dimensional] = myConvertUnit(reconSurface)
    % Kiểm tra giá trị lớn nhất để xác định đơn vị
%     if (max(reconSurface(:))-min(reconSurface(:))) < 10
    if max(abs(reconSurface(:))) > 1000
        scaleFactor = 10^3; % Chuyển từ nanomet sang micromet
        dimensional = 'micromet';
    else
        scaleFactor = 1;  
        dimensional = 'nanomet';
    end

    % Chuyển đổi đơn vị
    reconSurface = reconSurface / scaleFactor;
end

