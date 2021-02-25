%% Caclulating Second Reflections (Only works for one Tx).
Rx.SecondRefRSSI = zeros(size(Rx.xyz,1),1);
for i = 1:size(Rx.xyz,1)
    Rx.SeconReflWallJRSSI = zeros(size(wall.xyz1,1),1);
    for j = 1:size(wall.xyz1,1)
        %Initializing Parameters for Tx to First Reflection point.logs the wall between Tx and first reflection point on wall J
        Tx2FirstReflPintIntersecWalls = zeros(size(wall.xyz1,1),1); 
        % logs the Trans coeff of the wall between the Tx and the first reflection point
        Tx2FirstReflPintTransCoeff = ones(size(wall.xyz1,1),1); 
        %Initializing Parameters for First to Second Reflection ponit
        first2SecondReflPintIntersecWalls = zeros(size(wall.xyz1,1),1);
        first2SecondReflPintTransCoeff = ones(size(wall.xyz1,1),1);
        %Initializing Parameters for SECOND to Rx path
        second2RxIntersecWalls = zeros(size(wall.xyz1,1),1);
        second2RxReflPintTransCoeff = ones(size(wall.xyz1,1),1);
        Rx.SecondReflWallKRSSI = zeros(size(wall.xyz1,1),1);
        for k = 1:size(Tx.secondReflecWallj.xyz,1)
            if (sum(Tx.secondReflecWallj.xyz(k,:,j) ~= Tx.xyz) ~= 0) % checks if the Tx.secondReflecWallj lies on the Tx
                TxSecondRef2Rx.vec.xyz = Rx.xyz(i,:) - Tx.secondReflecWallj.xyz(k,:,j);
                TxSecondRef2Rx.dist = sqrt(sum(TxSecondRef2Rx.vec.xyz.^2,2));
                % Find intersection of Tx.secondReflecWallj with wall k
                 TxSecndRef2wallKIntd = dot(wall.xyz1(k,:) - Tx.secondReflecWallj.xyz(k,:,j), wall.normal.xyz(k,:),2) ./...
                     dot(TxSecondRef2Rx.vec.xyz, wall.normal.xyz(k,:),2);        
                 %check if there is intersection between the TxSecondreflection and the wall K
                 if (TxSecndRef2wallKIntd < 1 && TxSecndRef2wallKIntd > 0) 
                     secondReflectPointK = TxSecndRef2wallKIntd .* TxSecondRef2Rx.vec.xyz + Tx.secondReflecWallj.xyz(k,:,j);
                     % now check if secondReflectPointK actually lies on the finite plane(wall) K
                     if (prod(wall.minMax.x(k,:) - secondReflectPointK(1,1),2) < eps) && (prod(wall.minMax.y(k,:) -...
                             secondReflectPointK(1,2),2) < eps) && (prod(wall.minMax.z(k,:) - secondReflectPointK(1,3),2) < eps)
                         % At this point there is a path for second reflection, now chekc if there is a valid
                        % path for first reflection (LOS between first reflection and secondReflectPointK intersects with wall j
                         TxRefj2secondReflectPointK.vec.xyz = secondReflectPointK - Tx.wallReflec.xyz(j,:);
                         % Find the intersection of wall j with TxRef2secondReflectPointK
                         TxRefj2SecondReflectPointKWalljIntd = dot(wall.xyz1(j,:) - Tx.wallReflec.xyz(j,:),...
                             wall.normal.xyz(j,:),2) ./ dot(TxRefj2secondReflectPointK.vec.xyz, wall.normal.xyz(j,:),2); 
                         % check if there is intersection
                         if (TxRefj2SecondReflectPointKWalljIntd < 1 && TxRefj2SecondReflectPointKWalljIntd > 0) 
                             % now that there is intersection
                             firstReflecPointj = TxRefj2SecondReflectPointKWalljIntd .* TxRefj2secondReflectPointK.vec.xyz +...
                                 Tx.wallReflec.xyz(j,:);
                             % check if the intersection lies on a finite plane
                             if (prod(wall.minMax.x(j,:) - firstReflecPointj(1,1),2) < eps) && (prod(wall.minMax.y(j,:) -...
                                     firstReflecPointj(1,2),2) < eps)&&(prod(wall.minMax.z(j,:)-firstReflecPointj(1,3),2)<eps) 
                                 % At this point there is a path for second & first reflections,
                                 % 1- Find the reflection coefficient for both first and second reflections
                                 % 2- Count the walls between reflection paths
                                 % 1- Finding Reflection Coefficient
                                 tempSecondReflecAngle = acosd(abs(dot(TxSecondRef2Rx.vec.xyz,wall.normal.xyz(k,:),2) ./...
                                     ((sqrt(sum(TxSecondRef2Rx.vec.xyz.^2,2)) .* sqrt(sum(wall.normal.xyz(k,:).^2))))));
                                 tempFirstReflecAngle  = acosd(abs(dot(TxRefj2secondReflectPointK.vec.xyz,...
                                     wall.normal.xyz(j,:),2) ./ ((sqrt(sum(TxRefj2secondReflectPointK.vec.xyz.^2,2)) .*...
                                     sqrt(sum(wall.normal.xyz(j,:).^2))))));
                                 % Second Reflection factors baised on wall K for second reflections
                                if  k < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % if panel is a wall
                                    tempSecondReflecCoeff = wall.TE.refFac(k,(round(tempSecondReflecAngle)+1));
                                else % if panel is either ceiling or floor
                                    tempSecondReflecCoeff = wall.TM.refFac(k,(round(tempSecondReflecAngle)+1));
                                end
                                % First Reflection factors baised on wall j for second reflections
                                 if  j < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % if panel is a wall
                                    tempFirstReflecCoeff = wall.TE.refFac(j,(round(tempFirstReflecAngle)+1));
                                 else % if panel is either ceiling or floor
                                    tempFirstReflecCoeff = wall.TM.refFac(j,(round(tempFirstReflecAngle)+1));
                                 end
% The path of second reflection breaks down into 3 parts. Tx to First Reflection point. First to second
            % Reflection point. Second reflection point to RX. For each path, walls in between need ot be checked
                                 Tx2FirstReflPoint.vec.xyz = firstReflecPointj - Tx.xyz;
                                 firstReflPoint2SecondReflPoint.vec.xyz  = secondReflectPointK - firstReflecPointj;
                                 secondReflPoint2Rx.vec.xyz = Rx.xyz(i,:) - secondReflectPointK;
                                 % checking number of walls between TX and firstReflectionPoint
                                 for l = 1:size(wall.xyz1,1) % checking number of walls between TX and firstReflectionPoint
                                     if (l ~= j) % unecessary, only for safety measures as intersection D for same wall is zero.
                                         wallLIntdTx2FirstReflPoint = dot(wall.xyz1(l,:) - Tx.xyz, wall.normal.xyz(l,:),2) ./...
                                             dot(Tx2FirstReflPoint.vec.xyz, wall.normal.xyz(l,:),2);         
                                         if (wallLIntdTx2FirstReflPoint < 1 && wallLIntdTx2FirstReflPoint > 0)
                                             wallLintTx2Tx2FirstReflPoint = wallLIntdTx2FirstReflPoint .*...
                                                 Tx2FirstReflPoint.vec.xyz + Tx.xyz;
                                             if (prod(wall.minMax.x(l,:) - wallLintTx2Tx2FirstReflPoint(1,1),2) < eps) &&...
                                                     (prod(wall.minMax.y(l,:) - wallLintTx2Tx2FirstReflPoint(1,2),2)<eps)&&...
                                                     (prod(wall.minMax.z(l,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                              % now wall L is in between. Find the angle of incidence and transmission coefficient
                                                 Tx2FirstReflPintIntersecWalls(l,1) = 1; % logging the wall in between
                                                 tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),...
                                                     Tx2FirstReflPoint.vec.xyz,2)./(sqrt(sum(wall.normal.xyz(l,:).^2,2)).*...
                                                     sqrt(sum(Tx2FirstReflPoint.vec.xyz.^2,2))))); % finds the angle between
                                                if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                                                    % Transmission Coeffs for the intercepting walls between TX and refl point
                                                    Tx2FirstReflPintTransCoeff(l) = wall.TE.transFac(round...
                                                            (tempWallInterceptingAngle(l))+1);
                                                else % if panel is either ceiling or floor
                                                    Tx2FirstReflPintTransCoeff(l) = wall.TM.transFac(round(...
                                                            tempWallInterceptingAngle(l))+1);
                                                end 
                                             end 
                                         end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).
                                     end % if (l ~= j)
                                 end % for l = 1:size(wall.xyz1,1)
                                 % checking number of walls between the FIRSTReflectionPoint and SECONDreflection point
                                 for l = 1:size(wall.xyz1,1) % First Reflection Pint to Second
                                     wallLIntdFirstRefl2SecondReflPoint = dot(wall.xyz1(l,:) - firstReflecPointj,...
                                         wall.normal.xyz(l,:),2) ./ dot(firstReflPoint2SecondReflPoint.vec.xyz,...
                                         wall.normal.xyz(l,:),2);         
                                     if (wallLIntdFirstRefl2SecondReflPoint < 1 && wallLIntdFirstRefl2SecondReflPoint > 0)
                                         wallLIntFirstRefl2SecondReflPoint = wallLIntdFirstRefl2SecondReflPoint .*...
                                             firstReflPoint2SecondReflPoint.vec.xyz + firstReflecPointj;
                                         if (prod(wall.minMax.x(l,:) - wallLIntFirstRefl2SecondReflPoint(1,1),2) < eps)...
                                                 && (prod(wall.minMax.y(l,:) - wallLIntFirstRefl2SecondReflPoint(1,2),2)<eps)...
                                                 && (prod(wall.minMax.z(l,:) - wallLIntFirstRefl2SecondReflPoint(1,3),2) < eps)
                        % now wall L is in between. Find the angle of incidence and transmission coefficient
                                             first2SecondReflPintIntersecWalls(l,1) = 1; % logging the wall in between
                                             tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),...
                                                 firstReflPoint2SecondReflPoint.vec.xyz,2)./...
                                            (sqrt(sum(wall.normal.xyz(l,:).^2,2)) .* sqrt(sum...
                                            (firstReflPoint2SecondReflPoint.vec.xyz.^2,2))))); % finds the angle between
                                        % Transmission Coeffs for the intercepting walls between TX and refl point
                                            if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                                                first2SecondReflPintTransCoeff(l) = wall.TE.transFac(round...
                                                        (tempWallInterceptingAngle(l))+1);
                                            else % if panel is either ceiling or floor
                                                first2SecondReflPintTransCoeff(l) = wall.TM.transFac(round...
                                                        (tempWallInterceptingAngle(l))+1);
                                            end 
                                         end 
                                     end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).                             
                                 end % for l = 1:size(wall.xyz1,1)
                                 % checking number of walls between the SECOND reflection point to the RX.
                                 for l = 1:size(wall.xyz1,1) % SECOND to RX
                                     % To Avoid the Same wall that the second Reflection is bouncing of! Although this is 
                                     % for safety really
                                     if (l ~= k) 
                                         wallLIntdSecondRef2Rx = dot(wall.xyz1(l,:) - Rx.xyz(i,:), wall.normal.xyz(l,:),2) ./...
                                             dot(secondReflPoint2Rx.vec.xyz, wall.normal.xyz(l,:),2);         
                                         if (wallLIntdSecondRef2Rx < 1 && wallLIntdSecondRef2Rx > 0)
                                             wallLIntSecondRef2Rx = wallLIntdSecondRef2Rx .* secondReflPoint2Rx.vec.xyz +...
                                                 secondReflectPointK;
                                             if (prod(wall.minMax.x(l,:) - wallLIntSecondRef2Rx(1,1),2) < eps) &&...
                                                     (prod(wall.minMax.y(l,:) - wallLIntSecondRef2Rx(1,2),2) < eps) && ...
                                                     (prod(wall.minMax.z(l,:) - wallLIntSecondRef2Rx(1,3),2) < eps)
                         % now wall L is in between. Find the angle of incidence and transmission coefficient
                                                 second2RxIntersecWalls(l,1) = 1; % logging the wall in between
                                                 tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),...
                                                     secondReflPoint2Rx.vec.xyz,2)./...
                                                (sqrt(sum(wall.normal.xyz(l,:).^2,2)) .* sqrt(sum...
                                                (secondReflPoint2Rx.vec.xyz.^2,2))))); % finds the angle between
                                            % Transmission Coeffs for the intercepting walls between TX and refl point
                                                if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                                                    second2RxReflPintTransCoeff(l) = wall.TE.transFac(round...
                                                            (tempWallInterceptingAngle(l))+1);
                                                else % if panel is either ceiling or floor
                                                    second2RxReflPintTransCoeff(l) = wall.TM.transFac(round...
                                                            (tempWallInterceptingAngle(l))+1);
                                                end 
                                             end 
                                         end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).
                                     end % if (l ~= k)
                                 end % for l = 1:size(wall.xyz1,1)
% There is reflection so find the antenna gain and beam departure angle, departure angle for reflections, is the angle between
% reflection point of the reflecting wall and the TX image.
% Elevation angle (between beam and Z plane not it's normal) -90<ele<90 degrees
                                depBeamAngle.Ele = asind(Tx2FirstReflPoint.vec.xyz(1,3) ./ sqrt(sum...
                                    (Tx2FirstReflPoint.vec.xyz.^2,2)));  
                                if isnan(depBeamAngle.Ele)
                                    depBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                end
                                depBeamAngle.Azi = atan2(Tx2FirstReflPoint.vec.xyz(1,2),Tx2FirstReflPoint.vec.xyz(1,1)) *...
                                    (180/pi); % Azimuth angle (between x and beam) -180<azi<180
                                if isnan(depBeamAngle.Azi)
                                    depBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                end
                                % between 1 to Resolution
                                depBeamAngle.AziIndex = ((depBeamAngle.Azi + 180)/360).* (antennaGainRes - 1) + 1; 
                                depBeamAngle.Zen = abs(90-depBeamAngle.Ele);  
                                % between 1 to antennaGainRes
                                depBeamAngle.ZenIndex = (depBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; 
                                % Also Calculating the Angle of Arrival
                                arrBeamAngle.Ele = asin(-secondReflPoint2Rx.vec.xyz(1,3) ./ sqrt(sum...
                                    (secondReflPoint2Rx.vec.xyz.^2,2)));
                                arrBeamAngle.Azi = atan2(-secondReflPoint2Rx.vec.xyz(1,2),-secondReflPoint2Rx.vec.xyz(1,1))...
                                    * (180/pi);
                                if isnan(arrBeamAngle.Ele)
                                    arrBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                end
                                % Zenith angle is calculated and used to find the antenna gain
                                arrBeamAngle.Zen = abs(90-arrBeamAngle.Ele);  
                                if isnan(arrBeamAngle.Azi)
                                    arrBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                end
                                % between 1 to antennaGainRes
                                arrBeamAngle.ZenIndex = (arrBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; 
                                % between 1 to Resolution
                                arrBeamAngle.AziIndex = ((arrBeamAngle.Azi + 180)/360).* (antennaGainRes - 1) + 1; 
                                % Calculating the Second Reflection of Wall 1:K for First Reflection being of wall J 
                                Rx.SecondReflWallKRSSI(k) = 10^((Tx.power - (FPSLRefLoss + 20*log10(4*pi*(TxSecondRef2Rx.dist)...
                                    .* freq ./ lightVel)) + (10*log10(prod(Tx2FirstReflPintTransCoeff)))...
                                + (10*log10(tempFirstReflecCoeff)) + (10*log10(prod(first2SecondReflPintTransCoeff)))...
                                + (10*log10(tempSecondReflecCoeff)) + (10*log10(prod(second2RxReflPintTransCoeff))) +...
                                (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex))) + ... 
                                (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex))))/10) ...
                                .* complex(cos(2*pi*freq*TxSecondRef2Rx.dist./lightVel) , sin(2*pi*freq*TxSecondRef2Rx.dist./...
                                lightVel));
                             end 
                         end 
                     end 
                 end %(TxSecndRef2wallKIntd < 1 && TxSecndRef2wallKIntd > 0)
                % check the validity of the first projection being on the wall
            end % if Tx.secondReflecWallj ~= Tx.xyz 
        end % for k = size(Tx.secondReflecWallj.xyz,1)
        Rx.SeconReflWallJRSSI(j) = sum(Rx.SecondReflWallKRSSI);
    end % for j = 1:size(Tx.secondReflecWallj.xyz,3)
    Rx.SecondRefRSSI(i,1) = sum(Rx.SeconReflWallJRSSI); 
end % for i = 1:size(Rx.xyz,1)