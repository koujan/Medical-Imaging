% questions 2.1 and 2.2
clear all;clc;  % pixel size or pixel spacing
M_R=dicomread('MAMMOGRAPHY_RAW');
M_R_info=dicominfo('MAMMOGRAPHY_RAW');
disp('MAMMOGRAPHY_RAW image detailes:');
size=M_R_info.Height*M_R_info.Width
disp('PixelSpacing');
M_R_info.PixelSpacing
disp('PatientName');
M_R_info.PatientName
disp('PatientBirthDate');
M_R_info.PatientBirthDate
disp('OtherPatientName');
M_R_info.OtherPatientName
disp('PatientAge');
M_R_info.PatientAge
disp('PatientSex');
M_R_info.PatientSex

%imtool(M_R);

M_P=dicomread('MAMMOGRAPHY_PRESENTATION');
M_P_info=dicominfo('MAMMOGRAPHY_PRESENTATION');
disp('MAMMOGRAPHY_PRESENTATION image detailes:');
size=M_P_info.Height*M_P_info.Width
disp('PixelSpacing');
M_P_info.PixelSpacing
disp('PatientName');
M_P_info.PatientName
disp('PatientBirthDate');
M_P_info.PatientBirthDate
disp('OtherPatientName');
M_P_info.OtherPatientName
disp('PatientAge');
M_P_info.PatientAge
disp('PatientSex');
M_P_info.PatientSex

%imtool(M_P);

MRI1=dicomread('MRI01');
MRI_info1=dicominfo('MRI01');
disp('MRI image detailes:');
size=MRI_info1.Height*MRI_info1.Width
disp('PixelSpacing');
MRI_info1.PixelSpacing
disp('PatientName');
MRI_info1.PatientName
disp('PatientBirthDate');
MRI_info1.PatientBirthDate
disp('PatientAge');
MRI_info1.PatientAge
disp('PatientSex');
MRI_info1.PatientSex


US=dicomread('ultrasound');
US_info=dicominfo('ultrasound');
disp('ultrasound image detailes:');
size=US_info.Width*US_info.Height
disp('PatientName');
US_info.PatientName
disp('PatientBirthDate');
US_info.PatientBirthDate
disp('PatientAge');
US_info.PatientAge
disp('PatientSex');
US_info.PatientSex
%imtool(US);

%% question 2.3&2.4&2.5(please change the directory below to the one in your PC)
clear all;clc;

dirc=dir('C:\Users\Mohammad\Documents\MATLAB\VIBOT 2nd\MI\MRI\MRI'); % the directory in which the images stored in my computer
vol=uint16(zeros(512,512,22));
for i=3:length(dirc) % first two elements of the dirc variable are special characters
    vol(:,:,i-2)=dicomread(dirc(i).name);
end
figure;imhist(vol(:));
title('Histogram of the MRI volume');

%%%axial
figure,imshow(vol(:,:,10));
title('axial orientation (slice 10)');
figure,imshow(vol(:,:,11));
title('axial orientation (slice 11)');

%%%coronal
a=uint16(zeros(512,22));
b=uint16(zeros(512,22));
a(:,:)=vol(255,:,:);
b(:,:)=vol(254,:,:);
a=imresize(a,'Scale' ,[1,3/0.3125]);
b=imresize(b,'Scale' ,[1,3/0.3125]);  %use image resize
figure,imshow(a);
title('coronal orientation (slice 255)');
figure,imshow(b);
title('coronal orientation (slice 254)');

% sagital view

c=uint16(zeros(512,22));
d=uint16(zeros(512,22));
c(:,:)=vol(:,255,:);
d(:,:)=vol(:,254,:);
c=imresize(c,'Scale' ,[1,3/0.3125]);
d=imresize(d,'Scale' ,[1,3/0.3125] );  %use image resize
figure,imshow(c);
title('sagita orientation (slice 255)');
figure,imshow(d);
title('sagital orientation (slice 254)');

%Mammography
M_R=dicomread('MAMMOGRAPHY_RAW');
M_P=dicomread('MAMMOGRAPHY_PRESENTATION');
%close all;
subplot(1,2,1);
imshow(M_P);
title(' MAMMOGRAPHY RAW');
subplot(1,2,2);
imshow(M_R);
title('MAMMOGRAPHY PRESENTATION');











