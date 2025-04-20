clc; clear; 

% Đường dẫn đến thư mục chứa các file CSV
folderPaths = {
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data obj wave 5 29-3\Data obj wave 5 29-3',  % OBJ
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data Ref 30-3-3\Data Ref 30-3-3'             % REF
};

% Chỉ số cột và hàng cần lấy dữ liệu
columnRange = [2, 38];       % Tổng cộng 37 cột
rowRange = [105, 133];       % Tổng cộng 29 hàng
numRows = rowRange(2) - rowRange(1) + 1;
numCols = columnRange(2) - columnRange(1) + 1;



% Xử lý OBJ và REF
averageMatrix1 = processCSVData(folderPaths{1}, rowRange, columnRange, numRows, numCols);
averageMatrix2 = processCSVData(folderPaths{2}, rowRange, columnRange, numRows, numCols);

wavelength = 633;
averageMatrix1 = averageMatrix1 * wavelength;
averageMatrix2 = averageMatrix2 * wavelength;

%%
% Tạo trục X và Y
% Tạo trục X và Y thực sự theo đơn vị mm
col_real = linspace(-2.7, 2.7, 37);  % Trục X (chiều ngang, ứng với cột)
row_real = linspace(-2.1, 2.1, 29);  % Trục Y (chiều dọc, ứng với hàng)

% Tạo lưới tọa độ từ trục thực
[X, Y] = meshgrid(col_real, row_real);

% Vẽ bề mặt tái tạo
figure;
mesh(X, Y, averageMatrix1);
shading interp; colormap jet; colorbar;
title('Bề mặt tái tạo sóng vật');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Wavefront (nm)');


figure;
mesh(X, Y, averageMatrix2);
shading interp; colormap jet; colorbar;
title('Bề mặt tái tạo của sóng tham chiếu');
xlabel('X(mm)'); ylabel('Y(mm)'); zlabel('Nanomet');

% Trừ bề mặt sóng
% averageMatrix = averageMatrix1 - averageMatrix2;
averageMatrix = averageMatrix1;

% Hiển thị bề mặt sau khi trừ
figure;
mesh(X, Y, averageMatrix);
shading interp; colormap jet; colorbar;
title('Sai khác giữa 2 bề mặt');
xlabel('X(mm)'); ylabel('Y(mm)'); zlabel('Nanomet');


% Lưu kết quả tạm thời
 save('main_tai_tao_pha_bu.mat');
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
