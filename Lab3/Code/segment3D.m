% Lab3: 3D  segmentation section
%note: please put this .m file in the same directory as the '3DSegmentation.mat.mat'
clear all;clc;
load('3DSegmentation.mat.mat');
TP=nnz(segmentation&mask);
FN=nnz(~segmentation&mask);
FP=nnz(segmentation&~mask);
TN=nnz(~segmentation&~mask);
TPR=TP/(TP+FN);
FPR=FP/(FP+TN);
%%%% Calculating the statistical parameters %%%%%%%%
positive_mask=nnz(mask);
positive_seg=nnz(segmentation);
jacc=TP/(positive_mask+positive_seg-TP);
Dice=2*TP/(positive_seg+positive_mask);

%%%% Hausdorff distance%%%%
s=strel('disk',1,0);
seg_ero=imerode(segmentation,s);
seg_boundary=segmentation-seg_ero;
mask_erod=imerode(mask,s);
mask_boundary=mask-mask_erod;
s=0;
p=0;
for i=1:size(segmentation,3) 
    [row_mask,col_mask]=find(mask_boundary(:,:,i)>0);
    [row_seg,col_seg]=find(seg_boundary(:,:,i)>0);
    if(~isempty(row_seg))
        seg_ones_ind(:,s+1:s+length(row_seg))=[row_seg';col_seg';i*ones(1,length(row_seg))];
        s=size(seg_ones_ind,2);
    end
    if(~isempty(row_mask))
        mask_ones_ind(:,p+1:p+length(row_mask))=[row_mask';col_mask';i*ones(1,length(row_mask))];
        p=size(mask_ones_ind);
    end
                
end
dist_seg=zeros(size(seg_ones_ind,2),1);
for i=1:size(seg_ones_ind,2) % calculate the distance from each point in the segmented image to all points in the mask
    dist=sqrt((seg_ones_ind(1,i)-mask_ones_ind(1,:)).^2+(seg_ones_ind(2,i)-mask_ones_ind(2,:)).^2+(seg_ones_ind(3,i)-mask_ones_ind(3,:)).^2);
    dist_seg(i)=min(dist);
end
dist_seg_mask=max(dist_seg);

dist_mask=zeros(size(mask_ones_ind,2),1);
for i=1:size(mask_ones_ind,2) % calculate the distance from each point in the mask image to all points in the segmented
    dist=sqrt((mask_ones_ind(1,i)-seg_ones_ind(1,:)).^2+(mask_ones_ind(2,i)-seg_ones_ind(2,:)).^2+(mask_ones_ind(3,i)-seg_ones_ind(3,:)).^2);
    dist_mask(i)=min(dist);
end
dist_mask_sig=max(dist_mask);
Hausdorff=max(dist_seg_mask,dist_mask_sig); 


