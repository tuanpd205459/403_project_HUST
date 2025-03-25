%% start

%% Xóa dữ liệu trong workspace và đóng tất cả các figure
clear;
close all;
clc;

%% Additional file (thêm file ảnh cần thiết - nếu nằm ngoài thư mục source code)
filePath = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Ảnh mẫu';
addpath(filePath);




%% Thêm đường dẫn tới file chứa hệ số Zernike
folderPath = "C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data Obj wave 2\Data Obj wave 2";
columnIndex = 4;    % Lấy dữ liệu từ cột 4
rowRange = [1, 66]; % Lấy dữ liệu từ hàng 1 đến 66

meanValues = reconstructZernike.computeMeanFromCSV(folderPath, columnIndex, rowRange);
% % Hiển thị kết quả
% disp('Giá trị trung bình cộng của cột 4 từ hàng 1 đến 66 trong tất cả các file CSV:');
% disp(meanValues);


%% Reconstructing wavefront from Zernike polys
 zernike_coeffs = meanValues;
 reconstructZernike.zernike_reconstruction(zernike_coeffs, 2.2, 0.633, 1); 

%  reconstructZernike.zernike_reconstruction(zernike_coeffs,grid_size);
% reconstructZernike.zernike_reconstruction_non_circle(zernike_coeffs,grid_size,grid_size);

%% Biến toàn cục (tham số đầu vào các hàm)
DPD = 25;
he_so = 1;
poly_order = 3;
maxIntensity =100000 ;  % Cường độ tối đa để hiển thị ảnh

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

filePath = 'nRBC (41).bmp';

hologram = processing.loadHologram(inputManual, filePath, folder_path);

%% processFourier: thực hiện các phép biến đổi Fourier và chọn vùng bậc +1
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
wavelength = 633e-9;
% DPD=1; % Do phong dai cua he quang
reconSurface = (unwrapped_Phase .* wavelength .* he_so) / (4*pi);
reconSurface = reconSurface * 10^6;     

offSet = 10;      
reconSurface = reconSurface(offSet:end-offSet,offSet:end-offSet);  % cắt/chọn vùng để vẽ đồ thị

temp_r = reconSurface;
%% Chuyển đổi đơn vị (sang nanomet nếu cần)
[reconSurface, dimensional] = processing.postProcess.myConvertUnit(reconSurface);

%% Xác định chiều vân (ngang/dọc)
detectFringe = processing.postProcess.detectFringeSobel(reconSurface);
disp(detectFringe);
if strcmpi(detectFringe, 'vân ngang')
    reconSurface = rot90(reconSurface); % Xoay 90 độ theo chiều dương
end

%% Post Processing
imagesc(reconSurface);
hold on;    
title('Mặt phẳng pha'); % Đặt tiêu đề cho hình ảnh

% Vẽ đường thằng cắt ngang
positionLine = processing.postProcess.myDrawLine();  

crossLine = processing.postProcess.myCrossSection(reconSurface, positionLine);

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
enableDisplayCrossSection3D = true;     %bật tắt mặt cắt ngang ảnh 3D

visualization.display3D(surfaceParams, true);


%% end