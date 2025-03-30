function meanValues = computeMeanFromCSV(folderPath, columnIndex, rowRange)
    % Hàm đọc dữ liệu từ nhiều file CSV và tính trung bình cộng
    % folderPath  - Đường dẫn đến thư mục chứa các file CSV
    % columnIndex - Chỉ số cột cần lấy dữ liệu
    % rowRange    - Khoảng hàng cần lấy dữ liệu (vector [startRow, endRow])
    
    % Lấy danh sách tất cả các file CSV trong thư mục
    csvFiles = dir(fullfile(folderPath, '*.csv'));
    
    % Số lượng file CSV tìm thấy
    numFiles = length(csvFiles);
    
    % Kiểm tra nếu không có file CSV nào
    if numFiles == 0
        error('Không tìm thấy file CSV nào trong thư mục: %s', folderPath);
    end
    
    % Xác định số lượng hàng cần lấy
    numRows = rowRange(2) - rowRange(1) + 1;
    
    % Khởi tạo cell array để lưu dữ liệu từ các file CSV
    allData = cell(1, numFiles);
    
    % Đọc từng file CSV
    for i = 1:numFiles
        filename = fullfile(folderPath, csvFiles(i).name); % Tạo đường dẫn đầy đủ
        data = readmatrix(filename); % Đọc dữ liệu từ file CSV
        
        % Kiểm tra kích thước dữ liệu để tránh lỗi
        if size(data, 1) >= rowRange(2) && size(data, 2) >= columnIndex
            allData{i} = data(rowRange(1):rowRange(2), columnIndex); % Lấy dữ liệu
        else
            warning('File %s không đủ kích thước.', filename);
            allData{i} = NaN(numRows, 1); % Gán NaN nếu file không đủ dữ liệu
        end
    end
    
    % Chuyển allData từ cell array sang ma trận (mỗi cột là một file)
    dataMatrix = cell2mat(allData);
    
    % Tính trung bình cộng theo từng hàng (giá trị trung bình của cột được chọn)
    meanValues = mean(dataMatrix, 2, 'omitnan'); % Bỏ qua giá trị NaN khi tính trung bình
end
