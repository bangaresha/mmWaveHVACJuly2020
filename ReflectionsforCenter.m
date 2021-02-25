function [tempFirst2SecPointTransCoeff]  = ReflectionsforCenter(t,s,wall,reflectPointList,first2SecondPointjList)
    first2SecondPointjWallsd = dot(wall.xyz1(s,:) - reflectPointList(s,:), wall.normal.xyz(s,:),2) ./...
         dot(first2SecondPointjList, wall.normal.xyz(s,:),2);         
    if (first2SecondPointjWallsd < 1 && first2SecondPointjWallsd > 0 && abs(first2SecondPointjWallsd - 1) > eps &&...
        not(first2SecondPointjWallsd < eps))
        first2SecondPointjWallsxyz = first2SecondPointjWallsd .* first2SecondPointjList + reflectPointList(s,:);
        if((prod(wall.minMax.x(s,:)- first2SecondPointjWallsxyz(1,1),2)<eps)&&(prod(wall.minMax.y(s,:)-...
                first2SecondPointjWallsxyz(1,2),2)<eps)&&(prod(wall.minMax.z(s,:)-first2SecondPointjWallsxyz(1,3),2)<eps))
                first2SecondPointIntersecWall(s) = 1; % At this point wall s is in between
                intercepWallsIncAngle.first2SecondPoint(s) = acosd(abs(dot(wall.normal.xyz(s,:),first2SecondPointjList,...
                2)./(sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(first2SecondPointjList.^2,2)))));
                if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1)
                    tempFirst2SecPointTransCoeff(s)=wall.TE.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
                else % if panel is a wall
                    tempFirst2SecPointTransCoeff(s)=wall.TM.transFac(round(intercepWallsIncAngle.first2SecondPoint(s))+1);
                end
        end
    else
        tempFirst2SecPointTransCoeff(s) = 0;
    end
    if t > 1
        ReflectionsforCenter(t-1,s,wall,reflectPointList,first2SecondPointjList)
    end
end