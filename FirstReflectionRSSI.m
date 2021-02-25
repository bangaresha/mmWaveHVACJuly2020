function [tempReflecCoeff,first2SecondPointj, first2SecondPointjDist,depBeamAngle, arrBeamAngle, first2SecondIntd, reflectPoint,reflecjRssi]  = FirstReflectionRSSI(i,j,t,RxXYZ,wall,TxReflectXYZList,Tx1XYZ,antennaGainRes,ceillFloor,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,k)
% reflectPoint = 10000*ones(1,3);
reflectPoint = ones(1,3);
first2SecondPointj = zeros(2,3);                                                                    
first2SecondPointjDist = zeros(2,1);
tempReflecCoeff = zeros(2,1);
depBeamAngle.Ele = 0;
depBeamAngle.Azi = 0;
depBeamAngle.AziIndex = 0;
depBeamAngle.Zen = 0;
depBeamAngle.ZenIndex = 0;
arrBeamAngle.Ele = 0;
arrBeamAngle.Azi = 0;
arrBeamAngle.AziIndex = 0;
arrBeamAngle.Zen = 0;
arrBeamAngle.ZenIndex = 0;   
reflecjRssi = 0;
TxReflectXYZ = TxReflectXYZList(:,1:3);
first2SecondVecXYZ = RxXYZ - TxReflectXYZ(j,:); %vector bw TxImage and Rx
%first2SecondDist = sqrt(sum(first2SecondVecXYZ.^2,2));
first2SecondIntd = dot(wall.xyz1(j,:) - TxReflectXYZ(j,:),wall.normal.xyz(j,:),2) ./ ...
    dot(first2SecondVecXYZ, wall.normal.xyz(j,:),2);
if (first2SecondIntd < 1 && first2SecondIntd > 0) 
    reflectPoint = first2SecondIntd .* first2SecondVecXYZ + TxReflectXYZ(j,:);
    [tempReflecCoeff,first2SecondPointj, first2SecondPointjDist,depBeamAngle, arrBeamAngle,reflectPoint2]  = distanceBwPoints(i,j,t,first2SecondVecXYZ,wall,reflectPoint,ceillFloor,RxXYZ,Tx1XYZ,antennaGainRes,TxReflectXYZList,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,k);
    [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff1] = furtherReflections(j,wall,Tx1XYZ,ceillFloor,first2SecondPointj(1,:));
    [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff2] = furtherReflections(j,wall,reflectPoint2(2,:),ceillFloor,first2SecondPointj(2,:));
    if size(tempReflecCoeff,1) == 1
        reflecjRssi = RssiFunc(t,Tx,FPSLRefLoss,first2SecondPointjDist(1,:),first2SecondPointjDist(2,:),Tx2FirstReflPintTransCoeff1,Tx2FirstReflPintTransCoeff2,1,tempReflecCoeff(1,:),1,depBeamAngle,arrBeamAngle,lightVel,Rx2TxRefl.dist(i,1,j,1),freq,TxAntennaGainAE,RxAntennaGainAE);
    else
        reflecjRssi = RssiFunc(t,Tx,FPSLRefLoss,first2SecondPointjDist(1,:),first2SecondPointjDist(2,:),Tx2FirstReflPintTransCoeff1,Tx2FirstReflPintTransCoeff2,1,tempReflecCoeff(1,:),tempReflecCoeff(2,:),depBeamAngle,arrBeamAngle,lightVel,Rx2TxRefl.dist(i,1,j,1),freq,TxAntennaGainAE,RxAntennaGainAE);
    end
end
end 