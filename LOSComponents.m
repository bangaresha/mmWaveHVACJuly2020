%% Calculating LOS Component
for i = 1:size(Tx.xyz,1)
    % beam angle is measured in relation to the origin, unit vectors of [1,0,0],[0,1,0],[0,0,1]. This depends on the 
    % orientation of the TX. Elevation angle (between beam and Z plane note it's normal) -90<ele<90 degrees
    losBeamAngle.Tx.Ele(:,i) = asind(-RxTx.vec.xyz(:,3,i) ./sqrt(sum(RxTx.vec.xyz(:,:,i).^2,2))); 
    % if nan turns the beam angle to 0
    losBeamAngle.Tx.Ele(find(isnan(losBeamAngle.Tx.Ele(:,i)) == 1),i) = 0; 
    % Azimuth angle (between x and beam) -180<azi<180
    losBeamAngle.Tx.Azi(:,i) = atan2(-RxTx.vec.xyz(:,2,i),-RxTx.vec.xyz(:,1,i))* (180/pi); 
    % Azimuth angle is unlikely to be nan due to use of atan2
    losBeamAngle.Tx.Azi(find(isnan(losBeamAngle.Tx.Azi(:,i)) == 1),i) = 0; 
    % between 1 to Resolution
    losBeamAngle.Tx.AziIndex(:,i) = ((losBeamAngle.Tx.Azi(:,i) + 180)./360).*(antennaGainRes - 1) + 1; 
    % Zenith angle is calculated and used to find the antenna gain
    losBeamAngle.Tx.Zen(:,i) = abs(90-losBeamAngle.Tx.Ele(:,i));  
    % between 1 to Resolution
    losBeamAngle.Tx.ZenIndex(:,i) = (losBeamAngle.Tx.Zen(:,i)./180).*(antennaGainRes - 1) + 1; 
    losBeamAngle.Rx.Ele(:,i) = -1 .* losBeamAngle.Tx.Ele(:,i);
    losBeamAngle.Rx.Zen(:,i) = abs(90 - losBeamAngle.Rx.Ele(:,i));
    losBeamAngle.Rx.ZenIndex(:,i) = (losBeamAngle.Rx.Zen(:,i)./180).* (antennaGainRes - 1) + 1;
    losBeamAngle.Rx.Azi(:,i) = atan2(RxTx.vec.xyz(:,2,i),RxTx.vec.xyz(:,1,i)) *(180/pi);
    losBeamAngle.Rx.Azi(find(isnan(losBeamAngle.Rx.Azi(:,i)) == 1),i) = 0; 
    losBeamAngle.Rx.AziIndex(:,i) = ((losBeamAngle.Rx.Azi(:,i) + 180)./360).*(antennaGainRes - 1) + 1; 
end
Tx2RxWalljd = zeros(size(wall.xyz1,1),1);
Tx2RxWalljxyz = zeros(size(wall.xyz1,1),3);
Tx2RxVec = zeros(size(Tx.xyz,1),3);
Tx2RxIntersectingWalls = zeros(size(Tx2RxWalljd));
Rx.LosRssi = zeros(size(Rx.xyz,1),1);
% figure
TxRxRayList = [];
for k = 1:size(Rx.xyz,1)  
    Tx2RxVec(1,:) = Rx.xyz(k,:) - Tx.xyz(1,:); 
    rayX = [Tx.xyz(1,1) Rx.xyz(k,1)];
    rayY = [Tx.xyz(1,2) Rx.xyz(k,2)];
    rayZ = [Tx.xyz(1,3) Rx.xyz(k,3)];
    TxRxRayList = [TxRxRayList; rayX rayY rayZ];
    Tx2RxIntersectingWalls = zeros(size(Tx2RxWalljd));
    incidentAngle = zeros(size(Tx2RxWalljd));
    tempFresnelCoeff = ones(size(Tx2RxWalljd));
    for j = 1:size(wall.xyz1,1) 
        % find intersection with each wall and validate it . % Scalar value of the line between TX & Rx
        Tx2RxWalljd(j) = dot(wall.xyz1(j,:)-Tx.xyz(1,:),wall.normal.xyz(j,:),2)./dot(Tx2RxVec(1,:),wall.normal.xyz(j,:),2); 
        if (Tx2RxWalljd(j)<1 && Tx2RxWalljd(j)>0)
            % Intersection point with wall j
            Tx2RxWalljxyz(j,:) = Tx2RxWalljd(j) .* Tx2RxVec(1,:) + Tx.xyz(1,:); 
            if (prod(wall.minMax.x(j,:)-Tx2RxWalljxyz(j,1),2)<eps)&&(prod(wall.minMax.y(j,:)...
                -Tx2RxWalljxyz(j,2),2)<eps)&&(prod(wall.minMax.z(j,:)-Tx2RxWalljxyz(j,3),2)<eps)
                Tx2RxIntersectingWalls(j) = 1; % At this point the intersection is definite
                % Angle between the beam and intersecting wall
                incidentAngle(j) = acosd(abs(dot(wall.normal.xyz(j,:),Tx2RxVec(1,:),2)./...
                    (sqrt(sum(wall.normal.xyz(j,:).^2,2)) .*sqrt(sum(Tx2RxVec(1,:).^2,2)))));
                rayX = [Tx.xyz(1,1) Tx2RxWalljxyz(j,1) Rx.xyz(k,1)];
                rayY = [Tx.xyz(1,2) Tx2RxWalljxyz(j,2) Rx.xyz(k,2)];
                rayZ = [Tx.xyz(1,3) Tx2RxWalljxyz(j,3) Rx.xyz(k,3)];
                if j < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % checks if it's wall or ceiling
                    % Invoking Temporary Transmission Coefficients for walls
                    tempFresnelCoeff (j) = wall.TE.transFac(j,(round(incidentAngle(j))+1)); 
                else
                    % Invoking Temporary Transmission Coefficients for Ceiling
                    tempFresnelCoeff (j) = wall.TM.transFac(j,(round(incidentAngle(j))+1));  
                end
            end
        end
    end
%     fill3(wall.X, wall.Y, wall.Z,wall.C)
%     hold on
%     plot3(rayX,rayY,rayZ); %,'LineStyle','none','Marker','.','Color','Black');
    Rx.LosRssi(k)=Rx.LosRssi(k)+10.^((Tx.power(1)-(FPSLRefLoss+20.*log10(4*pi*((RxTx.dist(k,1,1)>=refDistance).*...
        RxTx.dist(k,1,1)+(RxTx.dist(k,1,1)<refDistance)*refDistance).* freq./lightVel)-10.*log10(prod(tempFresnelCoeff)))+...
       (TxAntennaGainAE(round(losBeamAngle.Tx.AziIndex(k,1)),round(losBeamAngle.Tx.ZenIndex(k,1))))+...
       (RxAntennaGainAE(round(losBeamAngle.Rx.AziIndex(k,1)),round(losBeamAngle.Rx.ZenIndex(k,1)))))/10) .*...
       complex(cos(2*pi*freq*RxTx.dist(k,1,1)./lightVel + pi) ,sin(2*pi*freq*RxTx.dist(k,1,1)./lightVel + pi));
end
% hold off
TxRxX = [];
TxRxY = [];
TxRxZ = [];
LOSCount = 0;
figure
f = fill3(wall.X, wall.Y, wall.Z,wall.C);
hold on
for l=1:length(TxRxRayList)
    LOSCount = LOSCount + 1;
    TxRxX = [TxRxX; TxRxRayList(l,1) TxRxRayList(l,2)];
    TxRxY = [TxRxY; TxRxRayList(l,3) TxRxRayList(l,4)];
    TxRxZ = [TxRxZ; TxRxRayList(l,5) TxRxRayList(l,6)];
    plot3([TxRxRayList(l,1) TxRxRayList(l,2)], [TxRxRayList(l,3) TxRxRayList(l,4)], [TxRxRayList(l,5) TxRxRayList(l,6)]);
end
title("Line of Sight Only");
hold off
alpha(f, 0.5);