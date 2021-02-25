RayTracingStructure

%% EQUATING THE PANELS (WALLS) IN 3D - Find Walls Normals
wall.normal.xyz = (cross(wall.xyz2 - wall.xyz1,wall.xyz3 - wall.xyz1,2));
wall.unitNormal.xyz = wall.normal.xyz ./ repmat(sqrt(sum(wall.normal.xyz.^2,2)),1,3);
Tx.Proj.xyz(:,:,1) = repmat((dot((wall.xyz1 - repmat(Tx.xyz(1,:),size(wall.xyz1,1),1)),wall.unitNormal.xyz,2)./dot(...
    wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz+ repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1);
Tx.Reflec.xyz(:,:,1) = repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1)+ 2.*...
    (Tx.Proj.xyz(:,:,1) - repmat(Tx.xyz(1,:),size(wall.unitNormal.xyz,1),1));
TxProjXYZList = Tx.Proj.xyz(:,:,1);
TxReflectXYZList = Tx.Reflec.xyz(:,:,1);
tempRef.xyz(:,:,1) = Tx.Reflec.xyz(:,:,1);
%% Calculating the Image of Tx across each wall (every Tx.WallReflec should be images across all walls).
for t = 2
    [TxSecondProjWallj, TxSecondReflecWallj] = ImageTx(tempRef,wall);
    tempRef.xyz = TxSecondReflecWallj.xyz;
    TxProjXYZList = [TxProjXYZList TxSecondProjWallj.xyz(:,:,1) TxSecondProjWallj.xyz(:,:,2) TxSecondProjWallj.xyz(:,:,3) TxSecondProjWallj.xyz(:,:,4)];
    TxReflectXYZList = [TxReflectXYZList TxSecondReflecWallj.xyz(:,:,1) TxSecondReflecWallj.xyz(:,:,2) TxSecondReflecWallj.xyz(:,:,3) TxSecondReflecWallj.xyz(:,:,4)];
end
 
%% Distance Of TX(s) From Every Mesh Node (RXi), Its vector and unit vector
RxTx.vec.xyz(:,1:3,1) = repmat(Tx.xyz(1,:),size(Rx.xyz,1),1) - Rx.xyz;
RxTx.dist(:,1,1) = sqrt(sum(RxTx.vec.xyz(:,1:3,1).^2,2));
RxTx.unitVec.xyz(:,1:3,1) = RxTx.vec.xyz(:,1:3,1) ./ repmat(RxTx.dist(:,1,1),1,3);

%% Line Vector Between TxReflection & Rx
for j = 1:size(Tx.Reflec.xyz,1)
    Rx2TxRefl.vec.xyz(:,1:3,j,1)=repmat(Tx.Reflec.xyz(j,:,1),size(Rx.xyz,1),1) - Rx.xyz; %4th dimension represents Tx
    Rx2TxRefl.dist(:,1,j,1) = sqrt(sum(Rx2TxRefl.vec.xyz(:,1:3,j,1).^2,2)); 
    Rx2TxRefl.unitVec.xyz(:,1:3,j,1) = Rx2TxRefl.vec.xyz(:,1:3,j,1) ./repmat(Rx2TxRefl.dist(:,1,j,1),1,3);
end

LOSComponents
timestimes = 0;
RayTracingReflections
PropagationMaps2

% bandwidth = 80E6;
bandwidth = 2.16E9;
N0 = 1E-4;
capacity = bandwidth*log2(1 + max(Rx.TotalRSSILayer)/(N0*bandwidth));