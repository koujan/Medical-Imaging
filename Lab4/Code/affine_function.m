function [e]=affine_function(par,scale,Imoving,Ifixed,mtype,ttype)
% This function affine_registration_image, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% I=affine_registration_image(parameters,scale,I1,I2,type);
%
% input,
%   parameters (in 2D) : Rigid vector of length 3 -> [translateX translateY rotate]
%                        or Affine vector of length 7 -> [translateX translateY  
%                                           rotate resizeX resizeY shearXY shearYX]
%
%   parameters (in 3D) : Rigid vector of length 6 : [translateX translateY translateZ
%                                           rotateX rotateY rotateZ]
%                       or Affine vector of length 15 : [translateX translateY translateZ,
%                             rotateX rotateY rotateZ resizeX resizeY resizeZ, 
%                             shearXY, shearXZ, shearYX, shearYZ, shearZX, shearZY]
%   
%   scale: Vector with Scaling of the input parameters with the same length
%               as the parameter vector.
%   I1: The 2D/3D image which is affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   mtype: Metric type: s: sum of squared differences.
%
% outputs,
%   I: An image volume with the registration error between I1 and I2
%
% example,
%
% Function is written by D.Kroon University of Twente (July 2008)
x=par.*scale;

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

I3=affine_transform_2d_double(double(Imoving),double(M),1); % 3 stands for cubic interpolation
m=size(Ifixed,1);
n=size(Ifixed,2);
% metric computation (Here is the similarity measure which is "squared differences")
switch mtype
    case 'sd' %squared differences
        e=sum((I3(:)-Ifixed(:)).^2)/numel(I3);
    case 'm'
        I3=im2uint8(I3);
        Ifixed=im2uint8(Ifixed);
        e=0;
        PA=imhist(I3);PA=PA/(m*n);
        PB=imhist(Ifixed);PB=PB/(m*n);
        Iout=zeros(256,256);
        for i=1:size(Ifixed,1)
            for j=1:size(Ifixed,2)
                Iout(Ifixed(i,j)+1,I3(i,j)+1)=Iout(Ifixed(i,j)+1,I3(i,j)+1)+1;
            end
        end
        Iout=Iout./(m*n);
        for i=1:size(Iout,1)
            for j=1:size(Iout,2)
                if((PA(i)*PB(j))~=0 && Iout(i,j)~=0)
                    e=e+Iout(i,j)*log2(Iout(i,j)/(PA(i)*PB(j)) );
                end
            end
        end
        e=-e; % to be maximized by the optimizer
    otherwise
        error('Unknown metric type');
end;


