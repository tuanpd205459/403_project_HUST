%% Tính toán đa thức Zernike
function radial = zernike_radial(r,n,m)
    % Functions required for use: elliptical_crop
%     hàm tính toán đa thức Zernike
%     Đầu vào
%         r:      bán kính
%         n:      bậc quang sai
%         m:      bậc phương vị
%     Đầu ra:
%         giá trị Đa thức Zernike

    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0 || m < 0
        error('n and m must both be positive in radial function')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end
    if n == m
        radial = r.^n;
    elseif n - m == 2
        radial = n*zernike_radial(r,n,n)-(n-1)*zernike_radial(r,n-2,n-2);
    else
        H3 = (-4*((m+4)-2)*((m+4)-3)) / ((n+(m+4)-2)*(n-(m+4)+4));
        H2 = (H3*(n+(m+4))*(n-(m+4)+2)) / (4*((m+4)-1))  +  ((m+4)-2);
        H1 = ((m+4)*((m+4)-1) / 2)  -  (m+4)*H2  +  (H3*(n+(m+4)+2)*(n-(m+4))) / (8);
        radial = H1*zernike_radial(r,n,m+4) + (H2+H3 ./ r.^2).*zernike_radial(r,n,m+2);
        
        % Fill in NaN values that may have resulted from DIV/0 in prior
        % line. Evaluate these points directly (non-recursively) as they
        % are scarce if present.
        
        if sum(sum(isnan(radial))) > 0
            [row, col] = find(isnan(radial));
            c=1;
            while c<=length(row)
                x = 0;
                for k = 0:(n-m)/2
                    ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                    x = x + ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                end
                radial(row(c),col(c)) = x;
                c=c+1;
            end
        end

    end

end