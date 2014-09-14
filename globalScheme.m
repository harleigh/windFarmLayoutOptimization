%{
   globalScheme:
      Assuming the wind farm is a square, we build a globally optimal
      turbine placement when the wind is from 12 degrees only due north.

   Input:
     farmWidth: width of the square farm
     Nwt: number of turbines to be on the site
     dMin: relative minimum safe distance between turbines
     Trad: Radius of a turbine (center to propeller-tip)
     alpha: wake spreading constant

   Output:
     turbineSite: column vector of Nwt turbine sites such that
         * turbineSite(1:Nwt) are the x positions of the turbines
         * turbineSite(Nwt+1:2*Nwt) are the y positions of the turbines

   Notes:
    *A scheme is globally optimal (for this case where the wind is only
    comming from the north at 12 m/s) if no turbine is in another turbine's
    wake (and no voilation of the relative safe distance constraint of
    course)
    *The wind farm is a square, oriented (where L is the farm width)
                   (0,L)    (L,L)
                     *-------*
                     |       |
                     |       |
                     *-------*
                   (0,0)    (0,L)
    *Algorithm: We begin by placing one turbine at (0,L) the next turbine
    is going to be
            -Downwind by amount (dy)
            -To the right by ammount (dx)
    Such that
            -The relative distance between the two turbines is exactlyy
             dMin
            -The downwind turbine is not in the wake of the upwind turbine.
             Indeed, the blade of the downwind turbine is just outside of
             the wake zone of the upwind turbine
    Now, we know dx:  dx = wakeRadius(1,2)+Trad
                         = 2Trad + alpha*dwD(1,2)
    And we demand that the relative distance is dMin.  Hence we can find dy
    by an application of Pythagoream's theorem

Example on what we mean by Row 1 Row 2 Row 3: A zigZag, Row 1 is a zig, Row
2 is a zag and row 3 is a zig :)
               *      *
                *    *  *
                 *  *    *
                  *       *
               r1   r2  r3

%}
function turbineSite = globalScheme(farmWidth, Nwt, dMin, Trad, alpha)

%Quadratic formula, solving for dy in dy^2+(dy*alpha +2Trad)^2=dMin^2
%dy is the minimum yDistance between two turbines such that the relative
%distance between the turbines is equal to dMin.
a=1+alpha^2; b=4*Trad*alpha; c=4*Trad^2-dMin^2;
dy = (-b+sqrt(b^2-4*a*c))/(2*a); %dy is such that dy^2+(dy*alpha +2Trad)^2=dMin^2
dy = ceil(dy); %yDistance between turbines
dx = ceil(alpha*dy+2*Trad); %xDistance between turbines such that downwind turbine is not in upwind turbine's wake

yPosRow1 = (farmWidth:-dy:0)';
yPosRow2 = yPosRow1(end-1:-1:1);
yPosRow3 = yPosRow1(2:end);

ny = floor(length(yPosRow1)); %total # of turbines possible along y direction


xPosRow1 = (0:dx:dx*(ny-1))'; %ny number of x positions starting at x=0

%The first turbine on the second row needs to be offset by an ammount to
%avoid relative distance vilation with the second to last turbine on row 1
offset = max(dx,dMin-dx);
xPosRow2 = zeros(ny-1,1);
for i = 1:ny-1
    xPosRow2(i) = xPosRow1(end) + offset + (i-1)*dx;
end

%Calculation of xPositions for turbines in row 3
xPosRow3 = zeros(ny-1,1);
for i = 1:ny-1
    xPosRow3(i) = xPosRow2(end) + offset + (i-1)*dx;
end

xPos = [xPosRow1;xPosRow2;xPosRow3];
yPos = [yPosRow1;yPosRow2;yPosRow3];
turbineSite = [xPos(1:Nwt);yPos(1:Nwt)];

end
