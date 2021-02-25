function [TxReflection2Rx,first2SecondVecXYZ, first2SecondDist, first2SecondIntd, reflectPoint2,reflecjRssiT]  = ReflectionRSSI(i,j,t,RxXYZ,Tx,wall,TxReflectXYZ2,antennaGainRes,Tx1XYZ,ceillFloor,TxReflectXYZList,FPSLRefLoss,lightVel,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE)
    reflectPoint2 = ones(2,3);
    TxReflection2Rx = [];
    reflecjRssi = zeros(size(TxReflectXYZ2,1),1);
    % 1+3*(t-1)    (4 + 3*4)-1         (16 + 3*4)-1
    for k = 1:size(TxReflectXYZ2,1)
        %temp1 = TxReflectXYZ(:,((1 + 3*(j-1)):3*j));
        first2SecondVecXYZ = RxXYZ - TxReflectXYZ2(k,:,j); %vector bw TxImage and Rx
        first2SecondDist = sqrt(sum(first2SecondVecXYZ.^2,2));
        first2SecondIntd = dot(wall.xyz1(k,:) - TxReflectXYZ2(k,:,j),wall.normal.xyz(k,:),2) ./ ...
        dot(first2SecondVecXYZ, wall.normal.xyz(k,:),2);
        if (first2SecondIntd < 1 && first2SecondIntd > 0) 
            reflectPoint = first2SecondIntd .* first2SecondVecXYZ + TxReflectXYZ2(k,:,j);
            [tempReflecCoeff,first2SecondPointj, first2SecondPointjDist,depBeamAngle, arrBeamAngle,reflectPoint2] = distanceBwPoints(i,j,t,first2SecondVecXYZ,wall,reflectPoint,ceillFloor,RxXYZ,Tx1XYZ,antennaGainRes,TxReflectXYZList,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,k);
            [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff1] = furtherReflections(j,wall,Tx1XYZ,ceillFloor,first2SecondPointj(1,:));
            [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff2] = furtherReflections(j,wall,reflectPoint2(1,:),ceillFloor,first2SecondPointj(2,:));
            [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff3] = furtherReflections(j,wall,RxXYZ,ceillFloor,first2SecondPointj(3,:));
            reflecjRssi(k) = RssiFunc(t,Tx,FPSLRefLoss,first2SecondDist,1,Tx2FirstReflPintTransCoeff1,Tx2FirstReflPintTransCoeff2,Tx2FirstReflPintTransCoeff3,tempReflecCoeff(1,:),tempReflecCoeff(2,:),depBeamAngle,arrBeamAngle,lightVel,first2SecondDist,freq,TxAntennaGainAE,RxAntennaGainAE);
            TxReflection2Rx = [TxReflection2Rx; Tx.xyz(1,1) Tx.xyz(1,2) Tx.xyz(1,3); reflectPoint2(1,1) reflectPoint2(1,2) reflectPoint2(1,3); reflectPoint2(2,1) reflectPoint2(2,2) reflectPoint2(2,3); RxXYZ(1,1) RxXYZ(1,2) RxXYZ(1,3)];
        end
    end
 reflecjRssiT = sum(reflecjRssi); 
end 