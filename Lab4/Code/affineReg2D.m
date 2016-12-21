function [ Iregistered, M] = affineReg2D( Imoving, Ifixed,num )
%Example of 2D affine registration
%   Robert Marti  (robert.marti@udg.edu)
%   Based on the files from  D.Kroon University of Twente 

% clean
clear all; close all; clc;

% Read two imges 
Imoving=im2double(imread('lenag3.png')); 
Ifixed=im2double(imread('lenag2.png'));

% Smooth both images for faster registration
ISmoving=imfilter(Imoving,fspecial('gaussian'));
ISfixed=imfilter(Ifixed,fspecial('gaussian'));

mtype = 'sd'; % metric type: sd: ssd, m: mutual information, e: entropy 
ttype = 'r'; % rigid registration, options: r: rigid, a: affine
num=1;   % number of resolutions
% Parameter scaling of the Translation and Rotation
switch ttype
    case 'r'
        scale=[50 50 50];
        % Set initial affine parameters
        x=[0 0 0];
    case 'a'
        scale=[1 1 50 50 50 1 1];
        % Set initial affine parameters
        x=[0  0   pi/16    1    1   0  0];
    otherwise
        error('Unknown transformation type');
end
tic
for i=1:num
    x=x./scale; 
    x(1:2)=x(1:2)*2;  % for taking into account the translation for multi-resolution case
    [x]=fminsearch(@(x)affine_function(x,scale,imresize(ISmoving,1/(2^(num-i))),imresize(ISfixed,1/(2^(num-i))),mtype,ttype),x,optimset('Display','iter','MaxIter',1000,'PlotFcns',@optimplotfval, 'TolFun', 1.000000e-06,'TolX',1.000000e-06, 'MaxFunEvals', 1000*length(x)));
     % Scale the translation, resize and rotation parameters to the real values
    x=x.*scale;
end
toc
switch ttype
    case 'r'
         M=[ cos(x(3)) sin(x(3)) x(1);
            -sin(x(3)) cos(x(3)) x(2);
               0 0 1];
    case 'a'
         RT=[cos(x(3)) sin(x(3)) 0;
             -sin(x(3)) cos(x(3)) 0;
             0  0  1];
        sc=[x(4) 0 0;
              0 x(5) 0;
             0 0 1];
        sh=[1 x(6) 0;
           x(7) 1 0;
           0 0 1];
         T=[1 0 x(1);
          0 1 x(2);
             0 0  1];
         M=T*sh*sc*RT;
    otherwise
        error('Unknown transformation type');
  
end

% Transform the image 
Icor=affine_transform_2d_double(double(Imoving),double(M),0); % 3 stands for cubic interpolation

% Show the registration results
figure,
    subplot(2,2,1), imshow(Ifixed);
    subplot(2,2,2), imshow(Imoving);
    subplot(2,2,3), imshow(Icor);
    subplot(2,2,4), imshow(abs(Ifixed-Icor));
end

