function [tempReflecCoeff,first2SecondPointj, first2SecondPointjDist,depBeamAngle, arrBeamAngle,reflectPoint]  = distanceBwPoints(i,j,t,first2SecondVecXYZ,wall,reflectPoint2,ceillFloor,RxXYZ,Tx1XYZ,antennaGainRes,TxReflectXYZList,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,k)
%first2SecondVecXYZList((1+3*(j-1)):3*j)
reflectPoint = zeros(2,3);
tempReflecCoeff = zeros(2,1);
first2SecondPointj = zeros(3,3);
first2SecondPointjDist = 0;
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
if t == 1
    if(prod(wall.minMax.x(j,:)-reflectPoint2(1,1),2) < eps)&&(prod(wall.minMax.y(j,:)-...
                reflectPoint2(1,2),2)<eps)&& (prod(wall.minMax.z(j,:) - reflectPoint2(1,3),2) < eps)
        tempReflecAngle1 = acosd(abs(dot(first2SecondVecXYZ,wall.normal.xyz(j,:),2)./...
            ((sqrt(sum(first2SecondVecXYZ.^2,2)) .*sqrt(sum(wall.normal.xyz(j,:).^2))))));
        if  j < (size(wall.xyz1,1) -size(ceillFloor.xyz1,1) + 1)
            tempReflecCoeff = wall.TE.refFac(j,(round(tempReflecAngle1)+1));
        else % if panel is either ceiling or floor
            tempReflecCoeff = wall.TM.refFac(j,(round(tempReflecAngle1)+1));
        end
        first2SecondPointj = [(reflectPoint2 - Tx1XYZ); (RxXYZ - reflectPoint2)];
        first2SecondPointjDist = [sqrt(sum((reflectPoint2 - Tx1XYZ).^2,2)); sqrt(sum((RxXYZ - reflectPoint2).^2,2))];
        [depBeamAngle, arrBeamAngle] = anglesBwPoints(first2SecondPointj(1,:), first2SecondPointj(2,:),antennaGainRes);
    end
end
if t == 2
    if(prod(wall.minMax.x(k,:)-reflectPoint2(1,1),2) < eps)&&(prod(wall.minMax.y(k,:)-...
                reflectPoint2(1,2),2)<eps)&& (prod(wall.minMax.z(k,:) - reflectPoint2(1,3),2) < eps)
        tempReflecAngle2 = acosd(abs(dot(first2SecondVecXYZ,wall.normal.xyz(k,:),2)./...
            ((sqrt(sum(first2SecondVecXYZ.^2,2)) .*sqrt(sum(wall.normal.xyz(k,:).^2))))));
        if  k < (size(wall.xyz1,1) -size(ceillFloor.xyz1,1) + 1)
            tempReflecCoeff2 = wall.TE.refFac(k,(round(tempReflecAngle2)+1));
        else 
            tempReflecCoeff2 = wall.TM.refFac(k,(round(tempReflecAngle2)+1));
        end
        [tempReflecCoeff1,first2SecondPointj1, first2SecondPointjDist1,depBeamAngle1, arrBeamAngle1, first2SecondIntd1, reflectPoint1,reflecjRssi1] = FirstReflectionRSSI(i,j,1,reflectPoint2,wall,TxReflectXYZList,Tx1XYZ,antennaGainRes,ceillFloor,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,k);
        first2SecondPointj = [(reflectPoint1 - Tx1XYZ); (reflectPoint2 - reflectPoint1); (RxXYZ - reflectPoint2)];
        first2SecondPointjDist = [sqrt(sum((reflectPoint1 - Tx1XYZ).^2,2)); sqrt(sum((reflectPoint2 - reflectPoint1).^2,2)); sqrt(sum((RxXYZ - reflectPoint2).^2,2))];
        [depBeamAngle, arrBeamAngle] = anglesBwPoints(first2SecondPointj(1,:), first2SecondPointj(2,:),antennaGainRes);
        tempReflecCoeff = [tempReflecCoeff1; tempReflecCoeff2];
        reflectPoint = [reflectPoint1; reflectPoint2];
    end
end
end

%                 elseif t == 3
%                 first2SecondPointj = RxXYZ.xyz(i,:) - reflectPoint(j,:);
%                 first2SecondPointjDist = sqrt(sum(first2SecondPointj.^2,2));
%                 first2SecondPointj = [first2SecondPointj; first2SecondPointj];
%                 first2SecondPointjDist = [first2SecondPointjDist; first2SecondPointjDist];
%             else
%                 if j == 1
%                     first2SecondPointj = reflectPoint(j,:) - Tx1XYZ;
%                 else
%                     first2SecondPointj = reflectPoint(j,:) - reflectPoint(j-1,:);
%                 end
%                 first2SecondPointj = [first2SecondPointj; first2SecondPointj];
%                 first2SecondPointjDist = sqrt(sum(first2SecondPointj.^2,2));
%                 first2SecondPointjDist = [first2SecondPointjDist; first2SecondPointjDist];
%                 distanceBwPoints(t-1,i,first2SecondVecXYZ,wall,reflectPoint,ceillFloor,RxXYZ,Tx1XYZ)
