% Lab3: 2D  segmentation
%note: please put this .m file in the same directory as the examples folder
%which contains the mammogram images
clear all;clc;
jacc=zeros(6,4);
Dice=zeros(6,4);
Hausdorff=zeros(6,4);
ArUndCur=zeros(6,4);
Th_index=zeros(1,4);
for c=1:6 % loop to repeat the procedure for each segmentation algorithm
    cur_path_alg=['examples\alg' num2str(c) '\']; % directory of the current algorithm folder
    alg=dir([cur_path_alg '*.tif']);
    manual=dir(['examples\manual\' '*.tif' ]);
    Th=1:5:255;
    TPR=zeros(length(alg),length(Th));
    FPR=zeros(length(alg),length(Th));
    for i=1:length(alg) % loop to iterate on every image in the algi folder
        GT=imread(['examples\manual\' manual(i).name]);
        seg=imread([cur_path_alg alg(i).name]);
        temp=seg; %temporary parameter
        figure; 
        for j=1:length(Th) % loop to build the ROC of every algorithm for one image
            % thresholding step
            temp(seg>=Th(j))=255; 
            temp(seg<Th(j))=0;
            % calculating TPR, FPR
            TP=nnz(temp&GT);
            FN=nnz(~temp&GT);
            FP=nnz(temp&~GT);
            TN=nnz(~temp&~GT);
            TPR(i,j)=TP/(TP+FN);
            FPR(i,j)=FP/(FP+TN);
            scatter(FPR(i,j),TPR(i,j));hold on;
        end
        title(['TPR vs FPR for algorithm' num2str(c) ' and the clinical case number ' num2str(i)]);
        ArUndCur(c,i)=trapz(flip(FPR(i,:)),flip(TPR(i,:)));
        dist=sqrt((TPR(i,:)-1).^2+FPR(i,:).^2);
        [~,Th_index(i)]=min(dist);
        scatter(FPR(i,Th_index(i)),TPR(i,Th_index(i)),'fill');
    end
    best_Th=mean(Th(Th_index(:))); % best threshold for all images of one algorithm
    %%%% Calculating the statistical parameters %%%%%%%%
    for i=1:4 % for each image of a given algorithm
        GT=imread(['examples\manual\' manual(i).name]);
        seg=imread([cur_path_alg alg(i).name]);
        temp=seg;
        temp(seg>=best_Th)=255; 
        temp(seg<best_Th)=0;
        jacc(c,i)=nnz(temp&GT)/(nnz(GT)+nnz(temp)-nnz(temp&GT));
        Dice(c,i)=2*nnz(temp&GT)/(nnz(temp)+nnz(GT));
        %%%% Hausdorff distance%%%%
        s=strel('disk',1,0);
        temp_ero=imerode(temp,s);
        temp_boundary=temp-temp_ero;
        %figure,imshow(temp_boundary);
        GT_erod=imerode(GT,s);
        GT_boundary=GT-GT_erod;
        %figure,imshow(GT_boundary);
        [row_GT,col_GT]=find(GT_boundary>0);
        [row_temp,col_temp]=find(temp_boundary>0);
        if(isempty(row_temp)) % if the segmented image has zero TPR
           Hausdorff(c,i)=sqrt(size(temp,1)^2+size(temp,2)^2); % highest possible distance in this image
           continue;
        end
        dist_temp=zeros(length(row_temp),1);
        for r=1:length(row_temp)  % loop to calculate the distance form each point in the segmented image to the ground truth image
           dist=sqrt((row_temp(r)-row_GT).^2+(col_temp(r)-col_GT).^2);
           dist_temp(r)=min(dist);
        end
        dist_temp_GT=max(dist_temp);

        dist_GT=zeros(length(row_GT),1);
        for r=1:length(row_GT)   % loop to calculate the distance form each point in the ground truth image to the segmented image
           dist=sqrt((row_GT(r)-row_temp).^2+(col_GT(r)-col_temp).^2);
           dist_GT(r)=min(dist);
        end
        dist_GT_temp=max(dist_GT);
        Hausdorff(c,i)=max(dist_temp_GT,dist_GT_temp);
    end
end



