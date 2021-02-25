%% Caclulating Second Reflections (Only works for one Tx).
%Rx.reflecjRssi = zeros(4,1);
Tx1XYZ = Tx.xyz(1,:);
TxReflection2RxList = [];
RSSI3 = [];
RSSI2 = [];
impRespFirstReflect = [];
FirstReflectCount = 0;
figure
fill3(wall.X, wall.Y, wall.Z,wall.C)
hold on
for t = 1:2
    for i = 1:size(Rx.xyz,1)
        for j = 1:size(wall.xyz1,1)
            RxXYZ = Rx.xyz(i,:);
            TxReflectXYZ2 = reshape(TxReflectXYZList(:,4:15),[4,3,4]);
            if t == 1
                [tempReflecCoeff, first2SecondPointj, first2SecondPointjDist,depBeamAngle, arrBeamAngle, first2SecondIntd, reflectPoint,reflecjRssi1]  = FirstReflectionRSSI(i,j,t,RxXYZ,wall,TxReflectXYZList,Tx1XYZ,antennaGainRes,ceillFloor,FPSLRefLoss,lightVel,Tx,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE,0);
                Rx.reflecjRssi(j,1) = reflecjRssi1;
                RSSI2 = [RSSI2 reflecjRssi1];
                TxReflection1Rx(:,:,i) = [Tx.xyz(1,1) reflectPoint(1,1) Rx.xyz(i,1); Tx.xyz(1,2) reflectPoint(1,2) Rx.xyz(i,2); Tx.xyz(1,3) reflectPoint(1,3) Rx.xyz(i,3)];
                plot3([Tx.xyz(1,1) reflectPoint(1,1) Rx.xyz(i,1)], [Tx.xyz(1,2) reflectPoint(1,2) Rx.xyz(i,2)], [Tx.xyz(1,3) reflectPoint(1,3) Rx.xyz(i,3)]);
                % For 1st reflection
                firstReflectDist = first2SecondPointjDist(1,1) + first2SecondPointjDist(2,1);
                timeFirstReflect = (first2SecondPointjDist(1,1) + first2SecondPointjDist(2,1))/vel;
                %impRespFirstReflect = [impRespFirstReflect; reflecjRssi1*exp(-1i*firstReflectDist*pi)  timeFirstReflect];
                impRespFirstReflect = [impRespFirstReflect; abs(reflecjRssi1)  timeFirstReflect];
            else
                [TxReflection2Rx, first2SecondVecXYZ, first2SecondDist, first2SecondIntd, reflectPoint,reflecjRssiT]  = ReflectionRSSI(i,j,t,RxXYZ,Tx,wall,TxReflectXYZ2,antennaGainRes,Tx1XYZ,ceillFloor,TxReflectXYZList,FPSLRefLoss,lightVel,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE);
                Rx.reflecjRssi(j,:) = reflecjRssiT;
                RSSI3 = [RSSI3 reflecjRssiT];
                TxReflection2RxList = [TxReflection2RxList; TxReflection2Rx];
            end
        end 
        if t == 1
            Rx.ReflecRssi(i,:) = sum(sum(Rx.reflecjRssi,1),2);
        else
            Rx.ReflecRssi(i,:) = sum(Rx.reflecjRssi);
        end
    end
    RxReflect(t,:) = Rx.ReflecRssi;
end
title("First Relfection Only");
hold off

TxReflection12RxX = [];
TxReflection12RxY = [];
TxReflection12RxZ = [];
for l = 1:(size(TxReflection2RxList,1)/4)
    TxReflection12RxX = [TxReflection12RxX; TxReflection2RxList(l,1) TxReflection2RxList(l+1,1) TxReflection2RxList(l+2,1) TxReflection2RxList(l+3,1)];
    TxReflection12RxY = [TxReflection12RxY; TxReflection2RxList(l,2) TxReflection2RxList(l+1,2) TxReflection2RxList(l+2,2) TxReflection2RxList(l+3,2)];
    TxReflection12RxZ = [TxReflection12RxZ; TxReflection2RxList(l,3) TxReflection2RxList(l+1,3) TxReflection2RxList(l+2,3) TxReflection2RxList(l+3,3)];
end
figure
fill3(wall.X, wall.Y, wall.Z,wall.C)
hold on
plot3(TxReflection12RxX, TxReflection12RxY, TxReflection12RxZ);
title("Second Relfection Only");
hold off
alpha(f,.5)

% figure
% f = fill3(wall.X, wall.Y, wall.Z,wall.C);
% hold on
% plot3(TxReflection1Rx13, TxReflection1Rx23, TxReflection1Rx33);
% plot3(TxReflection12RxX, TxReflection12RxY, TxReflection12RxZ);
% title("Both Relfections");
% hold off
% alpha(f,.5)

RSSI1 = sum(10*log10(abs(Rx.LosRssi)))/length(Rx.LosRssi);
RSSI = RSSI2 + RSSI3;
RSSIdB = (10*log10(abs(RSSI)))/2;
RSSIAver = (((sum(RSSIdB)+RSSI1)/length(RSSIdB))/2) + 5
RSSI2Avg = sum(10*log10(abs(RSSI2)))/length(RSSI2);

figure
stem((abs(impRespFirstReflect(:,2)) - timeLOS),impRespFirstReflect(:,1));

% RSSIdB = -10*log10(abs(RSSI)) - Tx.power(1,:);
%             first2SecondPointIntersecWall  = zeros(size(wall.xyz1,1),1);
%             tempFirst2SecPointTransCoeff = ones(size(wall.xyz1,1),1);
%             for s = 1:size(wall.xyz1,1) 
%                 if t == 1 %|| Rx.reflecjRssi(s,t) == 0
%                     first2SecondPointjWallsd = dot(wall.xyz1(s,:) - Tx.xyz(1,:),wall.normal.xyz(s,:),2)./...
%                             dot(first2SecondPointj(1,:),wall.normal.xyz(s,:),2);
%                     if (first2SecondPointjWallsd < 1 && first2SecondPointjWallsd > 0 && abs(first2SecondPointjWallsd - 1) > eps)
%                         first2SecondPointjWallsxyz = first2SecondPointjWallsd .*first2SecondPointj(1,:) + Tx.xyz(1,:);                              
%                         if((prod(wall.minMax.x(s,:)- first2SecondPointjWallsxyz(1,1),2)<eps)&&(prod(wall.minMax.y(s,:)-...
%                              first2SecondPointjWallsxyz(1,2),2)<eps)&&(prod(wall.minMax.z(s,:)-first2SecondPointjWallsxyz(1,3),2)<eps))
%                             first2SecondPointIntersecWall = 1; % At this point wall s is in between
%                             intercepWallsIncAngle.first2SecondPoint(s) = acosd(abs(dot(wall.normal.xyz(s,:),...
%                                 first2SecondPointj(2,:),2)./(sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(first2SecondPointj(2,:).^2,2)))));
%                             if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
%                                 tempFirst2SecPointTransCoeff(s)=wall.TE.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
%                             else % if panel is a wall
%                                 tempFirst2SecPointTransCoeff(s)=wall.TM.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
%                             end
%                         end
%                     end
%                 elseif t == 3 
%                     %reflectPointList(:,(1+4*(s-1)):4*s,:)
%                     first2SecondPointjWallsd = dot(wall.xyz1(s,:) - reflectPoint(s,:),wall.normal.xyz(s,:),2)./...
%                         dot(first2SecondPointj(s,:),wall.normal.xyz(s,:),2);
%                     if (first2SecondPointjWallsd < 1 && first2SecondPointjWallsd > 0 && abs(first2SecondPointjWallsd - 1)>eps&&not(first2SecondPointjWallsd < eps))
%                         first2SecondPointjWallsxyz = first2SecondPointjWallsd .* first2SecondPointj(s,:) + reflectPoint(s,:);
%                         if(prod(wall.minMax.x(s,:) - first2SecondPointjWallsxyz(1,1),2) < eps) && (prod(wall.minMax.y(s,:) -...
%                             first2SecondPointjWallsxyz(1,2),2)<eps)&&(prod(wall.minMax.z(s,:) - first2SecondPointjWallsxyz(1,3),2)<eps)
%                                 first2SecondPointIntersecWall(s) = 1; % At this point wall s in between (reflection point to Rx)
%                                 intercepWallsIncAngle.first2SecondPoint(s) = acosd(abs(dot(wall.normal.xyz(s,:),first2SecondPointj(s,:),2)./...
%                                     (sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(first2SecondPointj(s,:).^2,2)))));
%                                 if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
%                                     tempFirst2SecPointTransCoeff(s)=wall.TE.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
%                                 else % if panel is either ceiling or floor
%                                     tempFirst2SecPointTransCoeff(s)=wall.TM.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
%                                 end
%                         end
%                     end
%                     Rx.reflecjRssi(s,t) = 10^((Tx.power(1,:) - (FPSLRefLoss + 20*log10(4*pi*(first2SecondPointjDist(s)) .* ...
%                     freq ./ lightVel)) + (10*log10(prod(tempFirst2SecPointTransCoeff(s))))+ (10*log10(tempReflecCoeff(s)))+...
%                     (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex(s)))) + ... 
%                     (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex(s)))))/10) ...
%                     .* complex(cos(2*pi*freq* first2SecondPointjDist(s)./lightVel + pi) ,...
%                     sin(2*pi*freq*first2SecondPointjDist(s)./lightVel + pi));
%             end 
        %                 fill3(wall.X, wall.Y, wall.Z,wall.C)
        %                 hold on
        %                 plot3([Tx.xyz(1,1) reflectPoint(1,1) Rx.xyz(k,1)], [Tx.xyz(1,2) reflectPoint(1,2) Rx.xyz(k,2)], ...
        %                         [Tx.xyz(1,3) reflectPoint(1,3) Rx.xyz(k,3)]);
 
