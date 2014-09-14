%{
   Function Name: trueCostEvaluation

   Purpose:
     Given Nwt turbine positions on a square wind farm, oriented as
                   (0,L)    (L,L)
                     *-------*
                     |       |
                     |       |
                     *-------*
                   (0,0)    (0,L)
     along with the wind directions and speed of wind along those wind
     directions, calculate the anual expected power output (in kilowatt
     hours, kWh)

   Inputs:
     x:  a column vector of turbine positions The turbine positions are
     ordered as, x(1:Nwt) are the xPositions, x(Nwt+1:2*Nwt) are teh
     yPositions.

   Outputs:
     expectedPower: Anual expected power output of the wind farm
                    configuration, in kWh.
     U: the onset wind at each turbine, (positions given in as the input x).
        That is, U represents the actual wind that each turbine is receiving.
        Hence the information for turbine 1 (x,y,U) is
             x(1)        (the xPosition of turbine 1)
             x(Nwt+1)    (the yPosition of turbine 1)
             x(2*Nwt+1)  (the onset wind speed for turbine 1 in the first wind
                         direction)
             x((w+2)*Nwd +1) onset wind for turbine 1 in the wth wind
                             direction

   Notes:
     When we speak of 'cost' we mean the expected anual power output of a
     wind farm configuration

     This function is very similar to the cost function used in the
     UsserFun, but differs in that here we are calculating the onset wind,
     whereas in the userFun the onset wind (for each wind direction) is a
     decision variable and hence is not directly computed.  Similarly
     whereas in the cost in the userFun, we allow the possibility of
     computing the cost under a modification (where being in wakes
     peanalize) here cost modification is always off, and hence we get the
     true cost of the wind farm layout wrt wind directions.
%}
function [ expectedPower, U ] = trueCostEvaluation( x ) %x is a column vector

    global u0 alpha Ct Nwt Trad D;
    global windDirections Nwd;
    
    %variable representing whether we are performing the cost with respect
    %to a cost modifier; since in this function we want the true 
    mode = 0;
    
    %column k is the onset wind speed for the turbines in wind direction k
    U = zeros(Nwt,Nwd);

    XLOC = 1;  %enumeration: location of xPositions in turbineProfiles is col 1
    YLOC = 2;  %enumeration: location of yPositions in turbineProfiles is col 2
    
    for w=1:Nwd
        %rotate the turbines into the wind direction (clock wise)
        rotatedTurbPos = [ x(1:Nwt), x(Nwt+1:2*Nwt)]*getCWRotationMatrix(windDirections(w));
        turbineProfiles = [ rotatedTurbPos(:,XLOC), rotatedTurbPos(:,YLOC) ];

        %Now that the turbines are rotated into the wind direction, we index
        %them according to down-wind and cross wind.  Due to the chosen grid
        %orientation (square: (0,L) (L,L) (0,0) (L,O) upLeft upRight lowLeft
        %lowRight respectivly) we sort the turbines in Decending order wrt yPos
        %and then ascending order wrt xPosition
        indexedTurbines = sortrows(turbineProfiles,[-YLOC,XLOC]);

        %Critticaly important: dwD and cwD are w.r.t. the sorted turbines
        dwD = calcDownWindDist(indexedTurbines(:,YLOC));
        cwD = calcCrossWindDist(indexedTurbines(:,XLOC));
        wakeRadius = calcWakeRadius(dwD, Trad, alpha );
        inWake = calcInWake( cwD, wakeRadius, Trad );
        [Q,~] = calcVelocityDeficit( cwD, dwD, inWake, wakeRadius, Trad, mode);
        
        %Note that the onset wind in the wind direction w is indexed with
        %respect to the indexed-Turbines.        
        U(:,w) = calcOnsetWind( Q, dwD, u0, D, alpha, Ct );
    end
    %Postcondition: U is an Nwt by Nwd matrix, column w in U is the onset
    %wind (with respect to the indexed turbine postions) along wind
    %direction w
    
    %since we only consider wind at 12 m/s the power curve is 0.3*U^3
    expectedPower = sum(U(:).^3)*0.3;
end

%
% end costEvaluation
%