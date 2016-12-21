clear all;clc;
M_R=dicomread('MAMMOGRAPHY_RAW');
M_P=dicomread('MAMMOGRAPHY_PRESENTATION');
subplot(2,3,1);
imshow(M_R);title('Raw image');

M_R(M_R>=5200)=0;
Correction = 65535 * (double(M_R)/65535).^0.35;  % gamma transformation
Correction=65536-uint16(Correction);subplot(2,3,2);imshow(Correction);title('gamma transformation');
Correction=imsharpen(Correction,'Radius',5);subplot(2,3,3);imshow(Correction);title('sharpening');
Correction=adapthisteq(Correction,'NumTiles', [9 9],'clipLimit',0.5,'Distribution','uniform');
Correction=65536-Correction;subplot(2,3,4);imshow(Correction);title('adaptive histogram equalization');
for i=1:size(Correction,1)
    for j=1000:size(Correction,2)
        if(Correction(i,j)>1)
            Correction(i,j)=65536-Correction(i,j);
        end
    end
end
subplot(2,3,5);imshow(Correction);title('colors inversion inside breast');

H = fspecial('average',[5 5]);
Correction = imfilter(Correction,H);
subplot(2,3,6);
imshow(Correction);title('smoothing');%figure;imshow(M_P);break;

