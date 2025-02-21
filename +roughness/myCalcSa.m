%% Hàm tính Sa
function Sa = myCalcSa (inputMatrix, poly_order)
    [numRows, numCols] = size(inputMatrix);
    z_xy = processing.postProcess.myDifference(inputMatrix, poly_order);   % Tính sai lệch các giá trị
    abs_z = abs(z_xy(:)); 
    Sa = sum(abs_z(:)) / (numRows * numCols); 
    
end