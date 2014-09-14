%{
   Given a degree, return a clock-wise rotation matrix (in radians)

   Notes:
     * Rose Maps are given in terms of Degrees (i.e.: not radians)
     * For the orientation of a Rose Map, consider a clock on a wall, or a
       wrist watch:
           0 deg on rose map is 12o'clock on a clock  NORTH
           30 deg on rose map is one o'clock
           60 deg on rose map is two o'clock
           90 deg on rose map is three o'clock        EAST
           ...
           180 deg on rose map is six o'clock         SOUTH
           ...
           270 deg on rose map is nine o'clock        WEST
           ...
           330 deg on rose map is eleven o'clock

       Hence, when considering wind comming from 90 degrees on a Rose Map,
       we need to rotate the wind farm 90 degrees counter clockwise in
       order for the wind of 90 degrees to come from due-north

%}
function [ clockwiseRotationMatrix ] = getCWRotationMatrix( degree )

    radCos=cos(degree*pi/180);
    radSin=sin(degree*pi/180);
    clockwiseRotationMatrix = [radCos radSin; -radSin radCos ];

end

