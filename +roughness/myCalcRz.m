%% Hàm tính Rz
function output = myCalcRz(inputSuface, inputMeanLine)
    differenceValue = inputSuface - inputMeanLine;
    sorted_differenceValue = sort(differenceValue(:));
    
    top5_sum = sum(sorted_differenceValue(end-4:end));
    bottom5_sum = sum(sorted_differenceValue(1:5));
    output = (top5_sum - bottom5_sum) / 5;

end
