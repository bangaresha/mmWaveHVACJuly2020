function [depBeamAngle, arrBeamAngle] = anglesBwPoints(first2SecondPointj1,first2SecondPointj2,antennaGainRes)
depBeamAngle.Ele = asind(first2SecondPointj1(1,3)./sqrt(sum(first2SecondPointj1(1,:).^2,2)));  
if isnan(depBeamAngle.Ele)
    depBeamAngle.Ele = 0; % if nan turns the beam angle to 0
end
depBeamAngle.Azi = atan2(first2SecondPointj1(1,2),first2SecondPointj1(1,1))*(180/pi); 
if isnan(depBeamAngle.Azi)
    depBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
end
depBeamAngle.AziIndex = ((depBeamAngle.Azi + 180)/360).*(antennaGainRes - 1) + 1; % between 1 to Resolution
depBeamAngle.Zen = abs(90-depBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
depBeamAngle.ZenIndex = (depBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; % between 1 to antennaGainRes
%depBeamAngleList = [depBeamAngleList depBeamAngle];
arrBeamAngle.Ele = asin(-first2SecondPointj2(1,3) ./sqrt(sum(first2SecondPointj2(1,:).^2,2)));
if isnan(arrBeamAngle.Ele)
    arrBeamAngle.Ele = 0; % if nan turns the beam angle to 0
end
arrBeamAngle.Azi = atan2(-first2SecondPointj2(1,2),-first2SecondPointj2(1,1))*(180/pi);
if isnan(arrBeamAngle.Azi)
    arrBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
end
arrBeamAngle.AziIndex = ((arrBeamAngle.Azi + 180)/360).*(antennaGainRes - 1) + 1; % between 1 to Resolution
arrBeamAngle.Zen = abs(90-arrBeamAngle.Ele); %Zenith angle is calculated and used to find the antenna gain
arrBeamAngle.ZenIndex = (arrBeamAngle.Zen /180).* (antennaGainRes - 1) + 1; % between 1 to antennaGainRes
%arrBeamAngleList = [arrBeamAngleList arrBeamAngle];
