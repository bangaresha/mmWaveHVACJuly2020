function reflecjRssi = RssiFunc(t,Tx,FPSLRefLoss,first2SecondPointjDist1,first2SecondPointjDist2,tempFirst2SecPointTransCoeff1,tempFirst2SecPointTransCoeff2,tempFirst2SecPointTransCoeff3,tempReflecCoeff1,tempReflecCoeff2,depBeamAngle,arrBeamAngle,lightVel,Rx2TxRefl,freq,TxAntennaGainAE,RxAntennaGainAE)
if depBeamAngle.AziIndex == 0
    reflecjRssi = 10^((Tx.power(1,:) - (FPSLRefLoss + 20*log10(4*pi*(first2SecondPointjDist1 + first2SecondPointjDist2) .* ...
                freq ./ lightVel)) + (10*log10(prod(tempFirst2SecPointTransCoeff1)))+ ...
                (10*log10(prod(tempFirst2SecPointTransCoeff2))) + (10*log10(prod(tempFirst2SecPointTransCoeff3))) +...
                (10*log10(tempReflecCoeff1))+ (10*log10(tempReflecCoeff2))+ 10*log10(1/(4060*1)) + (t*10*log10(1/first2SecondPointjDist1))+ ...
                (t*10*log10(1/first2SecondPointjDist2)))/10).* complex(cos(2*pi*freq* Rx2TxRefl./lightVel + pi),...
                sin(2*pi*freq*Rx2TxRefl./lightVel + pi));
else
    reflecjRssi = 10^((Tx.power(1,:) - (FPSLRefLoss + 20*log10(4*pi*(first2SecondPointjDist1 + first2SecondPointjDist2) .* ...
                freq ./ lightVel)) + (10*log10(prod(tempFirst2SecPointTransCoeff1)))+...
                (10*log10(prod(tempFirst2SecPointTransCoeff2))) + (10*log10(prod(tempFirst2SecPointTransCoeff3))) +...
                (10*log10(tempReflecCoeff1))+ (10*log10(tempReflecCoeff2)) + 10*log10(1/(4060*1)) +...
                (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex))) + ... 
                (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex))) + ...
                (t*10*log10(1/first2SecondPointjDist1)) + (t*10*log10(1/first2SecondPointjDist2)))/10) ...
                .* complex(cos(2*pi*freq* Rx2TxRefl./lightVel + pi) ,...
                sin(2*pi*freq*Rx2TxRefl./lightVel + pi));
end
end