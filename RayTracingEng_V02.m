RayTracingStructure

%% EQUATING THE PANELS (WALLS) IN 3D - Find Walls Normals
wall.normal.xyz = (cross(wall.xyz2 - wall.xyz1,wall.xyz3 - wall.xyz1,2));
wall.unitNormal.xyz = wall.normal.xyz ./ repmat(sqrt(sum(wall.normal.xyz.^2,2)),1,3);
%Findin Projecn of Tx on eachpanel https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
Tx.wallProj.xyz(:,:,1) = repmat((dot((wall.xyz1 - repmat(Tx.xyz(1,:),size(wall.xyz1,1),1)),wall.unitNormal.xyz,2)./dot(...
    wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz+ repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1);
% Calculating the reflection (mirror) of Tx accross each panel
Tx.wallReflec.xyz(:,:,1) = repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1)...
    + 2.* (Tx.wallProj.xyz(:,:,1) - repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1));

%% Calculating the Image of Tx across each wall (every Tx.WallReflec should be images across all walls).
for i = 1:size(wall.xyz1,1) % only works for first Tx
%Tx.secondProjWallj.xyz(:,:,x) contains the projections of wall x-th Tx.wallReflec.xyz(x,:,1)
    Tx.secondProjWallj.xyz(:,:,i) = repmat((dot((wall.xyz1 - repmat(Tx.wallReflec.xyz(i,:,1),size(wall.xyz1,1),1)),...
        wall.unitNormal.xyz,2)./dot(wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz + ...
        repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1); % only works for 1 TX so (:,:,1)
    Tx.secondReflecWallj.xyz(:,:,i) = repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1) + 2.*...
        (Tx.secondProjWallj.xyz(:,:,i) -repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1));
end

%% Distance Of TX(s) From Every Mesh Node (RXi), Its vector and unit vector
RxTx.vec.xyz(:,1:3,1) = repmat(Tx.xyz(1,:),size(Rx.xyz,1),1) - Rx.xyz;
RxTx.dist(:,1,1) = sqrt(sum(RxTx.vec.xyz(:,1:3,1).^2,2));
RxTx.unitVec.xyz(:,1:3,1) = RxTx.vec.xyz(:,1:3,1) ./ repmat(RxTx.dist(:,1,1),1,3);

%% Line Vector Between TxReflection & Rx
for j = 1:size(Tx.wallReflec.xyz,1)
    Rx2TxRefl.vec.xyz(:,1:3,j,1)=repmat(Tx.wallReflec.xyz(j,:,1),size(Rx.xyz,1),1) - Rx.xyz; %4th dimension represents Tx
    % 3rd dimension represents Tximage across the wall which the reflection took place
    Rx2TxRefl.dist(:,1,j,1) = sqrt(sum(Rx2TxRefl.vec.xyz(:,1:3,j,1).^2,2)); 
    Rx2TxRefl.unitVec.xyz(:,1:3,j,1) = Rx2TxRefl.vec.xyz(:,1:3,j,1) ./repmat(Rx2TxRefl.dist(:,1,j,1),1,3);
end

LOSComponents
timestimes = 0;
FirstReflections
SecondReflections
PropagationMaps

% bandwidth = 80E6;
bandwidth = 2.16E9;
N0 = 1E-4;
capacity = bandwidth*log2(1 + max(Rx.TotalRSSILayer)/(N0*bandwidth));