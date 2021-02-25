%% Calculating Multipath & Reflection Components
figure
for k = 1:size(Rx.xyz,1)
    Rx.reflecjRssi = zeros(size(wall.xyz1,1),1);
    for j = 1:size(wall.xyz1,1) 
        temp1 = Tx.wallReflec.xyz(j,:,1);
        TxRef2Rx.vec.xyz = Rx.xyz(k,:) - Tx.wallReflec.xyz(j,:,1); %vector bw TxImage and Rx
        TxRef2RxRefwallIntd = dot(wall.xyz1(j,:) - Tx.wallReflec.xyz(j,:,1),...
            wall.normal.xyz(j,:),2) ./ dot(TxRef2Rx.vec.xyz, wall.normal.xyz(j,:),2);
        if (TxRef2RxRefwallIntd < 1 && TxRef2RxRefwallIntd > 0) % d checks the reflection possibility from TxImage j
            reflectPointj = TxRef2RxRefwallIntd .* TxRef2Rx.vec.xyz +Tx.wallReflec.xyz(j,:,1);
            % xyz check the reflection possibility from TxImage j
            if(prod(wall.minMax.x(j,:)-reflectPointj(1,1),2) < eps)&&(prod(wall.minMax.y(j,:)-reflectPointj(1,2),2)<eps)&&...
                    (prod(wall.minMax.z(j,:) - reflectPointj(1,3),2) < eps)
%At this point there is a path for reflection, 1-Find the reflection coefficient.2-Count the walls between reflection paths
                % 1- Finding Reflection Coefficient
                tempReflecAngle=acosd(abs(dot(TxRef2Rx.vec.xyz,wall.normal.xyz(j,:),2)./...
                ((sqrt(sum(TxRef2Rx.vec.xyz.^2,2)) .*sqrt(sum(wall.normal.xyz(j,:).^2))))));
                if  j < (size(wall.xyz1,1) -size(ceillFloor.xyz1,1) + 1)% if panel is a wall
                    tempReflecCoeff = wall.TE.refFac(j,(round(tempReflecAngle)+1));
                else % if panel is either ceiling or floor
                    tempReflecCoeff = wall.TM.refFac(j,(round(tempReflecAngle)+1));
                end
      % 2- now that there is reflection, lets find the walls between the reflection paths
                Tx2ReflectPointj = reflectPointj - Tx.xyz(1,:);
                reflectPointj2Rx = Rx.xyz(k,:) - reflectPointj;
                Tx2ReflectPointjDist = sqrt(sum(Tx2ReflectPointj.^2,2));
                reflectPointj2RxDist = sqrt(sum(reflectPointj2Rx.^2,2));
                Tx2ReflectPointIntersecWall  = zeros(size(wall.xyz1,1),1);
                reflectPointj2RxIntersecWall = zeros(size(wall.xyz1,1),1);
% There is reflection so find the antenna gain and beam departure angle.Departure angle for reflections is the angle between
% reflection point on the wall and the TX image. Elevation angle(between beam and Z plane not it's normal) -90<ele<90
                depBeamAngle.Ele = asind(Tx2ReflectPointj(1,3) ./sqrt(sum(Tx2ReflectPointj.^2,2)));  
                if isnan(depBeamAngle.Ele)
                    depBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                end
                % Azimuth angle (between x and beam) -180<azi<180
                depBeamAngle.Azi = atan2(Tx2ReflectPointj(1,2),Tx2ReflectPointj(1,1))*(180/pi); 
                if isnan(depBeamAngle.Azi)
                    depBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                end
                depBeamAngle.AziIndex = ((depBeamAngle.Azi + 180)/360).*(antennaGainRes - 1) + 1; % between 1 to Resolution
                depBeamAngle.Zen = abs(90-depBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
                depBeamAngle.ZenIndex = (depBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; % between 1 to antennaGainRes
                % Also Calculating the Angle of Arrival
                arrBeamAngle.Ele = asin(-reflectPointj2Rx(1,3) ./sqrt(sum(reflectPointj2Rx.^2,2)));
                if isnan(arrBeamAngle.Ele)
                    arrBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                end
                arrBeamAngle.Azi = atan2(-reflectPointj2Rx(1,2),-reflectPointj2Rx(1,1))*(180/pi);
                if isnan(arrBeamAngle.Azi)
                    arrBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                end
                arrBeamAngle.AziIndex = ((arrBeamAngle.Azi + 180)/360).*(antennaGainRes - 1) + 1; % between 1 to Resolution
                arrBeamAngle.Zen = abs(90-arrBeamAngle.Ele); %Zenith angle is calculated and used to find the antenna gain
                arrBeamAngle.ZenIndex = (arrBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; % between 1 to antennaGainRes
                % a- Find Walls Between TX to Reflection point 
                tempTx2ReflpointTransCoeff = ones(size(wall.xyz1,1),1);
                tempReflpoint2RxTransCoeff = ones(size(wall.xyz1,1),1);
                for s = 1:size(wall.xyz1,1) 
                    % Finding Scalar value of intersection lines    
                    Tx2ReflectPointjWallsd = dot(wall.xyz1(s,:) - Tx.xyz(1,:),wall.normal.xyz(s,:),2)./...
                        dot(Tx2ReflectPointj,wall.normal.xyz(s,:),2);
                    reflectPointj2Rxd = dot(wall.xyz1(s,:) - reflectPointj,wall.normal.xyz(s,:),2)./...
                        dot(reflectPointj2Rx,wall.normal.xyz(s,:),2);
                    % Checking for finite plane intersection
                    if (Tx2ReflectPointjWallsd < 1 && Tx2ReflectPointjWallsd > 0 && abs(Tx2ReflectPointjWallsd - 1) > eps)
                        Tx2ReflectPointjWallsxyz = Tx2ReflectPointjWallsd .*Tx2ReflectPointj + Tx.xyz(1,:);                              
                        if((prod(wall.minMax.x(s,:)- Tx2ReflectPointjWallsxyz(1,1),2)<eps)&&(prod(wall.minMax.y(s,:)-...
                        Tx2ReflectPointjWallsxyz(1,2),2)<eps)&&(prod(wall.minMax.z(s,:)-Tx2ReflectPointjWallsxyz(1,3),2)<eps))
                            Tx2ReflectPointIntersecWall(s) = 1; % At this point wall s is in between
                            intercepWallsIncAngle.Tx2ReflPoint(s) = acosd(abs(dot(wall.normal.xyz(s,:),...
                            Tx2ReflectPointj,2)./(sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(Tx2ReflectPointj.^2,2)))));
                            % Transmission Coeffs for the intercepting walls between TX and refl point
                          if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                              tempTx2ReflpointTransCoeff(s)=wall.TE.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                          else % if panel is a wall
                              tempTx2ReflpointTransCoeff(s)=wall.TM.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                          end
                        end
                    end
% b- Find finite walls Between Reflection point and RX. if below has complicated rule, as the reflectPointj2Rxd tend to 
% be smaller than matlab's epsilon and sometimes a little bit smaller than epsilon but bigger than epsm
                    if (reflectPointj2Rxd<1&&reflectPointj2Rxd>0&&abs(reflectPointj2Rxd - 1)>eps&&not(reflectPointj2Rxd<eps))
                        reflectPointj2Rxxyz = reflectPointj2Rxd .* reflectPointj2Rx + reflectPointj;
                        if(prod(wall.minMax.x(s,:) - reflectPointj2Rxxyz(1,1),2) < eps) && (prod(wall.minMax.y(s,:) -...
                                reflectPointj2Rxxyz(1,2),2)<eps)&&(prod(wall.minMax.z(s,:) - reflectPointj2Rxxyz(1,3),2)<eps)
                            reflectPointj2RxIntersecWall(s,1) = 1; % At this point wall s in between (reflection point to Rx)
                            intercepWallsIncAngle.ReflPoint2Rx(s) = acosd(abs(dot(wall.normal.xyz(s,:),reflectPointj2Rx,2)./...
                            (sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(reflectPointj2Rx.^2,2)))));
                        % Transmission Coeffs for the intercepting walls between TX and refl point
                            if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                               tempReflpoint2RxTransCoeff(s)=wall.TE.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                            else % if panel is either ceiling or floor
                               tempReflpoint2RxTransCoeff(s)=wall.TM.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                            end
                        end
                    end    
                end % end of for S
% found number of walls between reflection paths. calculate the received signal at Rx from reflection point on wall j
% This Considers distance from Tx image to Rx at once 
                Rx.reflecjRssi(j,1) = 10^((Tx.power(1,:) - (FPSLRefLoss + 20*log10(4*pi*(Tx2ReflectPointjDist+...
                    reflectPointj2RxDist) .* freq ./ lightVel)) + (10*log10(prod(tempTx2ReflpointTransCoeff)))...
                    + (10*log10(tempReflecCoeff)) + (10*log10(prod(tempReflpoint2RxTransCoeff))) +...
                    (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex))) + ... 
                    (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex))))/10) ...
                    .* complex(cos(2*pi*freq* Rx2TxRefl.dist(k,1,j,1)./lightVel + pi) ,...
                    sin(2*pi*freq*Rx2TxRefl.dist(k,1,j,1)./lightVel + pi));
                TxReflection1Rx(:,:,k) = [Tx.xyz(1,1) reflectPointj(1,1) Rx.xyz(k,1); Tx.xyz(1,2) reflectPointj(1,2) Rx.xyz(k,2); ...
                        Tx.xyz(1,3) reflectPointj(1,3) Rx.xyz(k,3)];
%                 fill3(wall.X, wall.Y, wall.Z,wall.C)
%         hold on
%         plot3([Tx.xyz(1,1) reflectPointj(1,1) Rx.xyz(k,1)], [Tx.xyz(1,2) reflectPointj(1,2) Rx.xyz(k,2)], ...
%                         [Tx.xyz(1,3) reflectPointj(1,3) Rx.xyz(k,3)]);
           end % end of if reflection exist
        end % end of d check for reflection
    end% end of for j
    Rx.ReflecRssi(k,1) = sum(sum(Rx.reflecjRssi,1),2);
end
hold off