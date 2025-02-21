# Thuật toán tái tạo H-L-GAM
#### 🔹 **Mục đích**
Đoạn code này thực hiện **xử lý ảnh H và tính toán độ nhám bề mặt 3D**, bao gồm:  
✅ Đọc ảnh H  
✅ Biến đổi F và tháo gỡ pha (phase unwrapping)  
✅ Chuyển đổi đơn vị, cắt vùng ảnh  
✅ Xác định chiều vân, lấy mặt cắt ngang  
✅ Tính toán độ nhám **Ra, Rz, Sa, Sq, Sz**  
✅ Hiển thị đồ thị 2D & 3D  

---

## ⚙ **Cách Sử Dụng**

### 1️⃣ **Chuẩn bị môi trường**
Thêm thư mục chứa dữ liệu
- Nếu ảnh nằm ngoài thư mục chứa code, cần thêm đường dẫn vào MATLAB bằng lệnh:
```matlab

filePath = 'C:\Users\admin\Máy tính\Lab thầy Tùng\Tài liệu a Tuân\Ảnh mẫu';
addpath(filePath);
```
- Cập nhật đường dẫn ảnh **filePath** hoặc thư mục **folder_path** chứa ảnh hologram  

### 2️⃣ **Thiết lập tham số đầu vào**
Các biến có thể điều chỉnh:  
```matlab
DPD = 25;               % Hệ số phóng đại
he_so = 1;              % Hệ số tính toán chiều cao
poly_order = 3;         % Bậc đa thức fitting đường trung bình
maxIntensity = 100000;  % Ngưỡng cường độ tối đa của ảnh
```

### 3️⃣ **Đọc ảnh hologram**
Có hai cách:  
✔ **Bằng tay**: Nhập đường dẫn file cụ thể (`filePath`)  
✔ **Tự động**: Lấy ảnh mới nhất từ `folder_path`  
```matlab
inputManual = 0;   % 1 = Chọn file thủ công, 0 = Tự động lấy ảnh mới nhất
folder_path = 'C:\Users\admin\...\Thi nghiem\6-12';
filePath = '8.bmp';   
hologram = processing.loadHologram(inputManual, filePath, folder_path);
```

### 4️⃣ **Xử lý ảnh**
Thực hiện biến đổi Fourier và tháo gỡ pha:  
```matlab
wrappedPhase = processing.processFourier(hologram);
methodGroup = 'poisson';  % Chọn thuật toán tháo gỡ pha
unwrapped_Phase = unwrapping.unwrapPhase(wrappedPhase, methodGroup);
```
✅ **Thay đổi phương pháp tháo gỡ pha**  
Có thể đổi sang `ls`, `tie`, `2dweight`…  
```matlab
methodGroup = 'ls';
methodType = 'dct';
unwrapped_Phase = unwrapping.unwrapPhase(wrappedPhase, methodGroup, methodType);
```

### 5️⃣ **Tái tạo bề mặt 3D**
✔ Cắt bỏ viền ảnh để tối ưu hiển thị:  
```matlab
offSet = 10;
reconSurface = reconSurface(offSet:end-offSet, offSet:end-offSet);
```

✔ Chuyển đổi đơn vị:
```matlab
[reconSurface, dimensional] = processing.postProcess.myConvertUnit(reconSurface);
```

### 6️⃣ **Xác định chiều vân & mặt cắt ngang**
✔ **Xác định chiều vân (ngang/dọc)**  
```matlab
detectFringe = processing.postProcess.detectFringeSobel(reconSurface);
if strcmpi(detectFringe, 'vân ngang')
    reconSurface = rot90(reconSurface);
end
```
✔ **Vẽ đường cắt ngang và tính trung bình**  
```matlab
positionLine = processing.postProcess.myDrawLine();
crossLine = processing.postProcess.myCrossSection(reconSurface, positionLine);
meanLine = processing.postProcess.myMeanLine(crossLine, poly_order);
```

### 7️⃣ **Tính toán độ nhám 2D & 3D**
✔ **Độ nhám 2D:**  
```matlab
[Ra, Ra_line] = roughness.myCalcRa(crossLine, meanLine, DPD);
Rz = roughness.myCalcRz(crossLine, meanLine);
```
✔ **Độ nhám 3D:**  
```matlab
Sz = roughness.myCalcSz(reconSurface, poly_order);
Sa = roughness.myCalcSa(reconSurface, poly_order);
Sq = roughness.myCalcSq(reconSurface, poly_order);
```

✔ **Lưu thông số**  
```matlab
surfaceParams = {reconSurface, Ra, Rz, Sa, Sq, Sz, Ra_line, positionLine, crossLine, meanLine, dimensional, DPD};
```

### 8️⃣ **Hiển thị kết quả**
✔ **Biểu đồ 2D**  
```matlab
visualization.display2D(surfaceParams);
```
✔ **Biểu đồ 3D (mặt cắt ngang tùy chọn)**  
```matlab
enableDisplayCrossSection3D = true;   
visualization.display3D(surfaceParams, enableDisplayCrossSection3D);
```


✅ **Cấu trúc thư mục cần có**  
📂 **Project Folder**  
┣ 📂 `processing/` *(Chứa các hàm xử lý ảnh)*  
┣ 📂 `roughness/` *(Chứa các hàm tính độ nhám)*  
┣ 📂 `visualization/` *(Chứa hàm vẽ đồ thị)*  
┣ 📄 `main.m` *(File chính chạy chương trình)*  

---
