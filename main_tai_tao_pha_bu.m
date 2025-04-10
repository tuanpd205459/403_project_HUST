clc; clear; close all;

% Đường dẫn đến thư mục chứa các file CSV
folderPaths = {
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data obj wave 5 29-3\Data obj wave 5 29-3',  % OBJ
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data Ref 30-3-3\Data Ref 30-3-3'             % REF
};

% Chỉ số cột và hàng cần lấy dữ liệu
columnRange = [2, 38];
rowRange = [105, 133];
numRows = rowRange(2) - rowRange(1) + 1;
numCols = columnRange(2) - columnRange(1) + 1;



% Xử lý OBJ và REF
averageMatrix1 = processCSVData(folderPaths{1}, rowRange, columnRange, numRows, numCols);
averageMatrix2 = processCSVData(folderPaths{2}, rowRange, columnRange, numRows, numCols);

% Vẽ bề mặt OBJ và REF
[X, Y] = meshgrid(1:numCols, 1:numRows);

figure;
surf(X, Y, averageMatrix1);
shading interp; colormap jet; colorbar;
title('Bề mặt tái tạo OBJ');
xlabel('X'); ylabel('Y'); zlabel('Z');

figure;
surf(X, Y, averageMatrix2);
shading interp; colormap jet; colorbar;
title('Bề mặt tái tạo REF');
xlabel('X'); ylabel('Y'); zlabel('Z');

% Trừ bề mặt sóng
% averageMatrix = averageMatrix1 - averageMatrix2;
averageMatrix = averageMatrix1;

% Hiển thị bề mặt sau khi trừ
figure;
surf(X, Y, averageMatrix);
shading interp; colormap jet; colorbar;
title('Bề mặt sau khi trừ');
xlabel('X'); ylabel('Y'); zlabel('Z');


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
