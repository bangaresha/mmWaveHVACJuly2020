alpha1.x = (1 - mesh_.xNodeNum .* imageRSSIScale)./(boundary(1,1) - boundary(1,2));
beta.x  = 1 - (boundary(1,1) .* alpha1.x);
alpha1.y = (1 - mesh_.yNodeNum .* imageRSSIScale)./(boundary(2,1) - boundary(2,2));
beta.y  = 1 - (boundary(2,1) .* alpha1.y);
imageBoundary(1,:) = round(alpha1.x .* boundary(1,:) + beta.x);
imageBoundary(2,:) = round(alpha1.y .* boundary(2,:) + beta.y);
imageWalls.x=(reshape(round((alpha1.x .*[wall.xyz1(:,1);wall.xyz2(:,1);wall.xyz3(:,1);...
    wall.xyz4(:,1)])+beta.x),size(wall.xyz1,1),4));
imageWalls.y = (reshape(round((alpha1.y .* [wall.xyz1(:,2); wall.xyz2(:,2);wall.xyz3(:,2);...
    wall.xyz4(:,2)]) + beta.y),size(wall.xyz1,1),4));
imageTx.xy(:,1) = round((alpha1.x .* Tx.xyz(:,1)) + beta.x);
imageTx.xy(:,2) = round((alpha1.y .* Tx.xyz(:,2)) + beta.y);
structImage = false(imageBoundary(2,2),imageBoundary(1,2));
for i = 1:size(Tx.xyz,1)
    structImage(imageTx.xy(i,1),imageTx.xy(i,2)) = 1;
end
for j = 1:size(wall.xyz1,1)
    for i = 1:3
        [wallC,wallR] = bresenham(imageWalls.x(j,i),imageWalls.y(j,i),imageWalls.x(j,i+1),imageWalls.y(j,i+1)); 
        for k = 1:numel(wallC)
            structImage(wallC(k),wallR(k)) = 1;
        end
    end
end
structImage = structImage(1:imageBoundary(1,2),1:imageBoundary(2,2));
structImage = imrotate(structImage,90);
figure
imshow(structImage)
title("Structure")

%% Defining a Finite Panel (wall)
for i = 1:size(wall.xyz1,1)
    wall.minMax.x(i,:) = [min([wall.xyz1(i,1),wall.xyz2(i,1),wall.xyz3(i,1),wall.xyz4(i,1)]),...
        max([wall.xyz1(i,1),wall.xyz2(i,1),wall.xyz3(i,1),wall.xyz4(i,1)])];
    wall.minMax.y(i,:) = [min([wall.xyz1(i,2),wall.xyz2(i,2),wall.xyz3(i,2),wall.xyz4(i,2)]),...
        max([wall.xyz1(i,2),wall.xyz2(i,2),wall.xyz3(i,2),wall.xyz4(i,2)])];
    wall.minMax.z(i,:) = [min([wall.xyz1(i,3),wall.xyz2(i,3),wall.xyz3(i,3),wall.xyz4(i,3)]),...
        max([wall.xyz1(i,3),wall.xyz2(i,3),wall.xyz3(i,3),wall.xyz4(i,3)])];
end

%% 3D Formation of the Structure
wall.X = [wall.xyz1(:,1)';wall.xyz2(:,1)';wall.xyz3(:,1)';wall.xyz4(:,1)'];
wall.Y = [wall.xyz1(:,2)';wall.xyz2(:,2)';wall.xyz3(:,2)';wall.xyz4(:,2)'];
wall.Z = [wall.xyz1(:,3)';wall.xyz2(:,3)';wall.xyz3(:,3)';wall.xyz4(:,3)'];
wall.C = zeros(size(wall.X)); 

%% Claculating Fresnel Coefficients for Walls
for i = 1:size(wall.xyz1,1)
    [wall.TE.refFac(i,:),wall.TE.transFac(i,:),wall.TM.refFac(i,:),wall.TM.transFac(i,:)]...
        = FresnelCoefficients(1,wall.relativePerm(i,1),0:90,0);
end

%% Meshing The Boundary Volume
if numel(linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum)) == 1
    try
        zplaneHeight = Rx.xyz(1,3);
    catch
        zplaneHeight = str2num(str2mat(inputdlg('Please assign the RX simulation height:','Heigh Assignment')));
    end
    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),linspace(boundary(2,1),boundary(2,2),...
        mesh_.yNodeNum),zplaneHeight);
else
    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),linspace(boundary(2,1),boundary(2,2),...
        mesh_.yNodeNum),linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum));
end
% Rx.xyz = [reshape(X,[],1,1),reshape(Y,[],1,1),reshape(Z,[],1,1)];
% Rx.xyz = [reshape(X(21:30,21:30),[],1,1),reshape(Y(21:30,21:30),[],1,1),reshape(Z(21:30,21:30),[],1,1)];
Rx.xyz = [reshape(X(31:40,31:40),[],1,1),reshape(Y(31:40,31:40),[],1,1),reshape(Z(31:40,31:40),[],1,1)];

figure
fill3(wall.X, wall.Y, wall.Z,wall.C)
hold on
plot3(Rx.xyz(:,1),Rx.xyz(:,2),Rx.xyz(:,3),'LineStyle','none','Marker','.','Color','Black');
hold on
for i = 1:size(Tx.xyz,1)
    plot3(Tx.xyz(i,1),Tx.xyz(i,2),Tx.xyz(i,3),'LineStyle','none','Marker','*','Color','Red');
end
hold off
