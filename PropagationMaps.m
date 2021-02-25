%% Line Of Sight Propagation Map Only
Rx.TotalRSSI = 10*log10(abs(Rx.LosRssi));
Rx.TotalRSSI(find(isinf(Rx.TotalRSSI) == 1)) = 0;
imageRSSI = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalRSSILayer(:,i) = (Rx.TotalRSSI(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*...
        (mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageRSSI = mat2gray(reshape(Rx.TotalRSSILayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageRSSIScaled = imresize(imrotate(imageRSSI,90),imageRSSIScale);
    imageRSSIScaledOverlayed = imoverlay(imageRSSIScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageRSSIScaledOverlayed));
    colorbarLabels = min((Rx.TotalRSSILayer(:,i))) + (0:5) .* ((max(Rx.TotalRSSILayer(:,i))-...
        min(Rx.TotalRSSILayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['LOS Only @ Z-Plane height of ',num2str(zplaneHeight(i),'%10.2f'),'; Tx Power = ',num2str(Tx.power')]);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end
    
%% Reflection Propagation Map Only
Rx.TotalRSSI = 10*log10(abs(Rx.ReflecRssi));
Rx.TotalRSSI(find(isinf(Rx.TotalRSSI) == 1)) = 0;
Rx.TotalRSSI(find((Rx.TotalRSSI) == 0)) = min(min(Rx.TotalRSSI));
imageRSSI = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalRSSILayer(:,i) = (Rx.TotalRSSI(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):...
        (i*(mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageRSSI = mat2gray(reshape(Rx.TotalRSSILayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageRSSIScaled = imresize(imrotate(imageRSSI,90),imageRSSIScale);
    imageRSSIScaledOverlayed = imoverlay(imageRSSIScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageRSSIScaledOverlayed));
    colorbarLabels = min((Rx.TotalRSSILayer(:,i))) + (0:5) .* ((max(Rx.TotalRSSILayer(:,i))-...
        min(Rx.TotalRSSILayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['First Reflections Only@ Z-Plane height of ',num2str(zplaneHeight(i),'%10.2f'),...
        '; Tx Power = ',num2str(Tx.power')]);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end

%% Second Reflection Propagation Map Only
Rx.TotalRSSI = 10*log10(abs(Rx.SecondRefRSSI));
Rx.TotalRSSI(find(isinf(Rx.TotalRSSI) == 1)) = 0;
Rx.TotalRSSI(find((Rx.TotalRSSI) == 0)) = min(min(Rx.TotalRSSI));
imageRSSI = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalRSSILayer(:,i) = (Rx.TotalRSSI(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*...
        (mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageRSSI = mat2gray(reshape(Rx.TotalRSSILayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageRSSIScaled = imresize(imrotate(imageRSSI,90),imageRSSIScale);
    % To Take care of the Blackage
    blackageMask = (imageRSSI == 1);
    blackageMaskScaled = imresize(imrotate(blackageMask,90),imageRSSIScale);
    blackageMaskScaled = blackageMaskScaled | structImage;
    imageRSSIScaledOverlayed = imoverlay(imageRSSIScaled,blackageMaskScaled,[0,0,0]);
    imshow(rgb2gray(imageRSSIScaledOverlayed));
    colorbarLabels = min((Rx.TotalRSSILayer(:,i))) + (0:5) .* ((max(Rx.TotalRSSILayer(:,i))-...
        min(Rx.TotalRSSILayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['Second Reflections Only@ Z-Plane height of ',num2str(zplaneHeight(i),'%10.2f'),'; Tx Power = ',...
        num2str(Tx.power')]);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end
    
%% Reflection & Line Of Signt Propagation Map
Rx.TotalRSSI = 10*log10(abs(Rx.LosRssi + (reflectExaggerationFac * Rx.ReflecRssi) + (reflectExaggerationFac *...
    Rx.SecondRefRSSI) ));
imageRSSI = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end

for i = 1:mesh_.zNodeNum
    Rx.TotalRSSILayer(:,i) = (Rx.TotalRSSI(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*(mesh_.xNodeNum.*...
        mesh_.yNodeNum))));
    imageRSSI = mat2gray(reshape(Rx.TotalRSSILayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageRSSIScaled = imresize(imrotate(imageRSSI,90),imageRSSIScale);
    imageRSSIScaledOverlayed = imoverlay(imageRSSIScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageRSSIScaledOverlayed));
    colorbarLabels = min((Rx.TotalRSSILayer(:,i))) + (0:5) .* ((max(Rx.TotalRSSILayer(:,i))-...
        min(Rx.TotalRSSILayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['Z-Plane height of ',num2str(zplaneHeight(i),'%10.2f'),'; LOS = ',num2str(losFlag),'; Reflec = ',...
        num2str(reflectionFlag),'; Tx Power = ',num2str(Tx.power')]);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end