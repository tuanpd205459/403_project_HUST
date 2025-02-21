%% Hàm tính Sq
function Sq = myCalcSq (inputMatrix, poly_order)
    [numRows, numCols] = size(inputMatrix);
    z_xy = processing.postProcess.myDifference(inputMatrix, poly_order);   % Tính sai lệch các giá trị
    abs_z = abs(z_xy(:)); 
    Sq = sqrt(1/(numRows*numCols) * sum(abs_z(:).^2));
    
end