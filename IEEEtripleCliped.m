function [output_image] = IEEEtripleCliped(input_img)
% Implementation of "Triple clipped histogram-based medical image 
% enhancement using spatial frequency. 
% IEEE Transactions on NanoBioscience, 20(3), 278-286." - 
% Date of Publication: 04 March 2021 -------------
% Implementation by Waqar Mirza
% Input: input_img % I used a high dynamic medical image in dicom format 
% Output: output_image % The enhanced output image of the whole technique 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_img = double(dicomread('Subject_1.dcm'));
input_img = double(input_img);
img=linspace(min(input_img(:)),max(input_img(:)),256);
img=uint8(arrayfun(@(x) find(abs(img(:)-x)==min(abs(img(:)-x))),input_img));
%% Tripple clipped Dynamic Histogram Equalization TCDHE 
% Compute the number of rows and columns of the image
[rows, cols] = size(img);
% Compute the number of bins in the histogram
num_bins = 256;
% Divide the image into three equal sub-images with equal number of pixels
% in every sub image
sub_img1 = img(1:rows/3, :);
sub_img2 = img(rows/3+1:2*rows/3, :);
sub_img3 = img(2*rows/3+1:end, :);
% Compute the histogram of each sub-image
hist1 = imhist(sub_img1, num_bins);
hist2 = imhist(sub_img2, num_bins);
hist3 = imhist(sub_img3, num_bins);
% Compute the mean of each sub-histogram
mean_hist1 = mean(hist1);
mean_hist2 = mean(hist2);
mean_hist3 = mean(hist3);
% Clip each sub-histogram using the mean
hist1_clipped = min(hist1, mean_hist1);
hist2_clipped = min(hist2, mean_hist2);
hist3_clipped = min(hist3, mean_hist3);
% Compute the cumulative distribution function (CDF) of each clipped sub-histogram
cdf1 = cumsum(hist1_clipped) / sum(hist1_clipped);
cdf2 = cumsum(hist2_clipped) / sum(hist2_clipped);
cdf3 = cumsum(hist3_clipped) / sum(hist3_clipped);
% Compute the probability density function (PDF) of each clipped sub-histogram
pdf1 = hist1_clipped / sum(hist1_clipped);
pdf2 = hist2_clipped / sum(hist2_clipped);
pdf3 = hist3_clipped / sum(hist3_clipped);
% Compute the equalized histogram for each sub-image
eq_hist1 = uint8(255 * cdf1(sub_img1+1));
eq_hist2 = uint8(255 * cdf2(sub_img2+1));
eq_hist3 = uint8(255 * cdf3(sub_img3+1));
% Concatenate the equalized histograms to form the TCDHE image
tcdhe_out = [eq_hist1; eq_hist2; eq_hist3];
%%
% Display the input and tripple clipped dynamic histogram equalization image
% figure;
% subplot(1, 2, 1);
% imshow(input_img,[]);
% title('Input Image');
% subplot(1, 2, 2);
% imshow(tcdhe_out,[]);
% title('TCDHE Image');
%% % Apply DWT on tcdhe_out
[LL_tcdhe_out, LH_tcdhe_out, HL_tcdhe_out, HH_tcdhe_out] = dwt2(tcdhe_out, 'haar');
%%
LL_tcdhe_out = padarray(LL_tcdhe_out,[1,0],'post');
LH_tcdhe_out = padarray(LH_tcdhe_out,[1,0],'post');
HL_tcdhe_out = padarray(HL_tcdhe_out,[1,0],'post');
HH_tcdhe_out = padarray(HH_tcdhe_out,[1,0],'post');
%%
% Displaying the images 
% subplot(2,2,1), imshow(LL_tcdhe_out, []), title('LL')
% subplot(2,2,2), imshow(LH_tcdhe_out, []), title('LH')
% subplot(2,2,3), imshow(HL_tcdhe_out, []), title('HL')
% subplot(2,2,4), imshow(HH_tcdhe_out, []), title('HH')
% Here, LL_tcdhe_out, LH_tcdhe_out, HL_tcdhe_out, and HH_tcdhe_out represent 
% the approximation, horizontal detail, vertical detail, and diagonal detail 
% coefficients of tcdhe_out, respectively.
%% Apply DWT on input_image
[LL_input_image, LH_input_image, HL_input_image, HH_input_image] = dwt2(input_img, 'haar');
% Similarly, LL_input_image, LH_input_image, HL_input_image, and HH_input_image 
% represent the corresponding coefficients of input_image
%% Apply SVD on LL sub-band of tcdhe_out
[U_tcdhe_out,S_tcdhe_out,V_tcdhe_out] = svd(LL_tcdhe_out);
%% Apply SVD on LL sub-band of input_image
[U_input_image,S_input_image,V_input_image] = svd(LL_input_image);
%% Max --> U_tcdhe_out, V_tcdhe_out
% Compute the max singular value for U_tcdhe_out
max_singular_value_U_tcdhe_out = max(max(abs(U_tcdhe_out)));
% Compute the max singular value for V_tcdhe_out
max_singular_value_V_tcdhe_out = max(max(abs(V_tcdhe_out)));
%% Max --> U_input_image, V_input_image
% Compute the max singular value for U_input_image
max_singular_value_U_input_image = max(max(abs(U_input_image)));
% Compute the max singular value for V_input_image
max_singular_value_V_input_image = max(max(abs(V_input_image)));
%% Improvement Factor Calculation 
% Add the maximum singular values for U_tcdhe_out and V_tcdhe_out
max_singular_values_tcdhe_out = max_singular_value_U_tcdhe_out + max_singular_value_V_tcdhe_out;
% Add the maximum singular values for U_input_image and V_input_image
max_singular_values_input_image = max_singular_value_U_input_image + max_singular_value_V_input_image;
% Divide the sum of maximum singular values for tcdhe_out and input_image
max_singular_values_ratio = max_singular_values_tcdhe_out / max_singular_values_input_image;
%%  Summation N Calculation 
% I'm taking Beta value as  0.52 , According to paper any value between 0.05
% and 0.95 can be used.
be_ta = 0.05;
sum_N_l = (be_ta*max_singular_values_ratio*S_input_image);
%%
sum_N_r = (1-be_ta)* (1/max_singular_values_ratio) *(S_tcdhe_out);
%% 
sum_N = (sum_N_l + sum_N_r);
%% LL_n Calculation
% Pad the matrix 'U_tcdhe_out' with an extra row and column of zeros
ll_n = U_tcdhe_out*sum_N*V_tcdhe_out';
%% Spatial Frequency on high frequency components of tcdhe_out 
% Initialize spatial frequency variables
freq_LH_tcdhe_out = zeros(size(LH_tcdhe_out));
freq_HL_tcdhe_out = zeros(size(HL_tcdhe_out));
freq_HH_tcdhe_out = zeros(size(HH_tcdhe_out));
% Compute the spatial frequency of each high frequency component
components = {'LH_tcdhe_out', 'HL_tcdhe_out', 'HH_tcdhe_out'};
for i = 1:length(components)
    % Compute the 2D FFT of the component
    F = fft2(eval(components{i}));
    % Shift the zero-frequency component to the center of the spectrum
    F = fftshift(F);
    % Compute the magnitude of the frequency domain representation
    mag = abs(F);
    % Compute the spatial frequency
    [N, M] = size(eval(components{i}));
    dx = 1/M;
    dy = 1/N;
    fx = (-M/2:M/2-1)*dx;
    fy = (-N/2:N/2-1)*dy;
    [Fx, Fy] = meshgrid(fx, fy);
    freq = sqrt(Fx.^2 + Fy.^2);
    % Store the spatial frequency in the corresponding variable
    if strcmp(components{i}, 'LH_tcdhe_out')
        freq_LH_tcdhe_out = freq;
    elseif strcmp(components{i}, 'HL_tcdhe_out')
        freq_HL_tcdhe_out = freq;
    elseif strcmp(components{i}, 'HH_tcdhe_out')
        freq_HH_tcdhe_out = freq;
    end
end
%% Spatial Frequency on high frequency components of Input Image
% Initialize spatial frequency variables
freq_LH_input_image = zeros(size(LH_input_image));
freq_HL_input_image = zeros(size(HL_input_image));
freq_HH_input_image = zeros(size(HH_input_image));
% Compute the spatial frequency of each high frequency component
components_input = {'LH_input_image', 'HL_input_image', 'HH_input_image'};
for k = 1:length(components_input)
    % Compute the 2D FFT of the component
    F_in = fft2(eval(components_input{k}));
    % Shift the zero-frequency component to the center of the spectrum
    F_in = fftshift(F_in);
    % Compute the magnitude of the frequency domain representation
    mag_in = abs(F_in);
    % Compute the spatial frequency
    [N_in, M_in] = size(eval(components_input{k}));
    dx_in = 1/M_in;
    dy_in = 1/N_in;
    fx_in = (-M_in/2:M_in/2-1)*dx_in;
    fy_in = (-N_in/2:N_in/2-1)*dy_in;
    [Fx_in, Fy_in] = meshgrid(fx_in, fy_in);
    freq_in = sqrt(Fx_in.^2 + Fy_in.^2);
    % Store the spatial frequency in the corresponding variable
    if strcmp(components{i}, 'LH_input_image')
        freq_LH_input_image = freq_in;
    elseif strcmp(components{i}, 'HL_input_image')
        freq_HL_input_image = freq_in;
    elseif strcmp(components{i}, 'HH_input_image')
        freq_HH_input_image = freq_in;
    end
end
%%  Normalized SF on the high frequncy band images input image
SN_lh_i = freq_LH_input_image ./ (freq_LH_input_image + freq_LH_tcdhe_out);
SN_hl_i = freq_HL_input_image ./ (freq_HL_input_image + freq_HL_tcdhe_out);
SN_hh_i = freq_HH_input_image ./ (freq_HH_input_image + freq_HH_tcdhe_out);
%%  Normalized SF on the high frequncy band images input image
SN_lh_o = freq_LH_tcdhe_out ./ (freq_LH_input_image + freq_LH_tcdhe_out);
SN_hl_o = freq_HL_tcdhe_out ./ (freq_HL_input_image + freq_HL_tcdhe_out);
SN_hh_o = freq_HH_tcdhe_out ./ (freq_HH_input_image + freq_HH_tcdhe_out);
%% Image Fusion
% (LH_n)
LH_n = SN_lh_i .* LH_input_image + SN_lh_o.*LH_tcdhe_out;
% (HL_n)
HL_n = SN_hl_i .* HL_input_image + SN_hl_o.*HL_tcdhe_out;
% (HH_n)
HH_n = SN_hh_i .* HH_input_image + SN_hh_o.*HH_tcdhe_out;
%% Inverse Wavelet Transform to form the enhanced Medical Image
output_image = idwt2(ll_n,LH_n,HL_n,HH_n,'haar');
% Normalize the image to the range [0, 1]
output_image = (output_image - min(output_image(:))) ./ (max(output_image(:)) - min(output_image(:)));
% Scale the image to the range [0, 255]
output_image = output_image * 255;
% Convert the image to uint8 format
output_image = uint8(output_image);
imshow(output_image,[]);
end