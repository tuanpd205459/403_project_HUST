%% start

clear;
clc;

%% Đường dẫn đến thư mục chứa file CSV
folderPaths = {
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data obj wave 5 29-3\Data obj wave 5 29-3',  ... OBJ
    'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Data Ref 30-3-3\Data Ref 30-3-3'             % REF
};

columnIndex = 4;        % Cột cần đọc
rowRange = [1, 66];     % Hàng từ 1 đến 66

% Đọc dữ liệu cột 4, từ hàng 1 đến 66, và tính trung bình
meanVector1 = processCSVDataColumn(folderPaths{1}, columnIndex, rowRange);  % Vật
meanVector2 = processCSVDataColumn(folderPaths{2}, columnIndex, rowRange);  % Tham chiếu
meanVector1(1:3,:) = 0;
meanVector2(1:3,:) = 0;

meanValues = meanVector1 - meanVector2;  % song vat - song tham chieu

surface_zernike1 = reconstructZernike.zernike_reconstruction(meanVector1, 1.83, 633, 'Bề mặt tái tạo đa thức Zernike (sóng vật)'); 

%%
surface_zernike2 = reconstructZernike.zernike_reconstruction(meanVector2, 1.83, 633, 'Bề mặt tái tạo đa thức Zernike (sóng tham chiếu)'); 

surface_zernike3 = reconstructZernike.zernike_reconstruction(meanValues, 1.83, 633, 'Bề mặt tái tạo đa thức Zernike- sau khi trừ' ); 

save('main1_zernike.mat');



%% Hàm đọc cột dữ liệu và tính trung bình
function meanVector = processCSVDataColumn(folderPath, columnIndex, rowRange)
    csvFiles = dir(fullfile(folderPath, '*.csv'));
    numFiles = length(csvFiles);

    if numFiles == 0
        error('Không tìm thấy file CSV nào trong thư mục: %s', folderPath);
    end

    numRows = rowRange(2) - rowRange(1) + 1;
    dataMatrix = NaN(numRows, numFiles);

    for i = 1:numFiles
        data = readmatrix(fullfile(folderPath, csvFiles(i).name));
        [rows, cols] = size(data);

        % Kiểm tra kích thước hợp lệ
        if rows >= rowRange(2) && cols >= columnIndex
            dataMatrix(:, i) = data(rowRange(1):rowRange(2), columnIndex);
        else
            warning('File %s không đủ kích thước. Bỏ qua.', csvFiles(i).name);
        end
    end

    % Tính trung bình theo từng dòng (qua các file)
    meanVector = mean(dataMatrix, 2, 'omitnan');
% meanVector =dataMatrix;
end
