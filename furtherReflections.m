function [wallLIntdTx2FirstReflPoint, wallLintTx2Tx2FirstReflPoint, tempWallInterceptingAngle, Tx2FirstReflPintTransCoeff] = furtherReflections(j,wall,Tx1XYZ,ceillFloor,first2SecondPointj3)
wallLIntdTx2FirstReflPoint = zeros(1,3);
wallLintTx2Tx2FirstReflPoint = zeros(1,3);
tempWallInterceptingAngle = 0;
Tx2FirstReflPintTransCoeff = ones(size(wall.xyz1,1),1);
for l = 1:size(wall.xyz1,1) % checking number of walls between TX and firstReflectionPoint
     if (l ~= j) % unecessary, only for safety measures as intersection D for same wall is zero.
         wallLIntdTx2FirstReflPoint = dot(wall.xyz1(l,:) - Tx1XYZ, wall.normal.xyz(l,:),2) ./...
             dot(first2SecondPointj3, wall.normal.xyz(l,:),2);         
         if (wallLIntdTx2FirstReflPoint < 1 && wallLIntdTx2FirstReflPoint > 0)
             wallLintTx2Tx2FirstReflPoint = wallLIntdTx2FirstReflPoint .* first2SecondPointj3 + Tx1XYZ;
             if (prod(wall.minMax.x(l,:) - wallLintTx2Tx2FirstReflPoint(1,1),2) < eps) &&...
                     (prod(wall.minMax.y(l,:) - wallLintTx2Tx2FirstReflPoint(1,2),2)<eps)&&...
                     (prod(wall.minMax.z(l,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                 Tx2FirstReflPintIntersecWalls(l,1) = 1; 
                 tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),...
                     first2SecondPointj3,2)./(sqrt(sum(wall.normal.xyz(l,:).^2,2)).*...
                     sqrt(sum(first2SecondPointj3.^2,2))))); % finds the angle between
                if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) 
                    Tx2FirstReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                else 
                    Tx2FirstReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                end 
             end 
         end 
     end 
end 
end
 