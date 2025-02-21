%% Hàm tính Sz
function Sz_avg = myCalcSz (inputMatrix, poly_order)
    % Khởi tạo mảng để lưu kết quả Sz cho từng hàng
    Sz_row = zeros(size(inputMatrix, 1), 1);
    
    for i = 1:size(inputMatrix, 1)
        row_vector = inputMatrix(i, :);
        
        % Tính giá trị trung bình của hàng i
        meanLine_row = processing.postProcess.myMeanLine(row_vector, poly_order);
        
        % Tính độ lệch bình phương
        Rz_row = roughness.myCalcRz(row_vector, meanLine_row);
        
        % Tính tổng và căn bậc hai của tổng cho hàng i
        sum_Rz_row = sum(Rz_row);
        Sz_row(i) = sqrt(sum_Rz_row);
    end
    
    % Tính giá trị trung bình của Sz cho các hàng
    Sz_avg = mean(Sz_row);
    
    disp(['Sz theo các hàng: ', num2str(Sz_avg)]);
end