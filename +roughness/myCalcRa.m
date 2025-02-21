%% Hàm tính Ra
function [Ra, Ra_line] = myCalcRa(inputLine, input_meanline, DPD)
    if isrow(inputLine)
        inputLine = inputLine'; 
    end
    if isrow(input_meanline)
        input_meanline = input_meanline';  
    end
    lr = length(inputLine) * 3.45/DPD;  % Length of the surface in micrometers
    u = linspace(1, lr, length(inputLine));  
    
    % tính Ra = tích phân hình thang
    Ra = (1/lr) * trapz(u, abs(inputLine - input_meanline));
    
    % Kết quả
  %  disp(['Surface Roughness Ra: ', num2str(Ra), ' micrometers']);
    
    Ra_line = Ra + input_meanline;
end