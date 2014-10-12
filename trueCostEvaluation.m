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
     U: the onset wind of each turbine for each wind direction and each
        slice.  U is an Nwt by Nwd*Nws matrix such that
             U = [ U_w=1 U_w=2 ... U_w=Nwd ]
        where U_w=i is the onset wind of the Nwt turbines along the ith
        wind direction (the underscore is to signify subscript of course)
        Now U_w=i is a Nwt by Nws marix
             U_w=i = [U^s=1 U^s=2   U^s=Nws]
        Where U^s=j is the onset wind wrt the jth wind slice of the free
        wind and is a row vector of Nwt by 1
             U^s=j = [U_1 U_2 ... U_Nwt]
         where U_i is the ith wind turbine

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

    global alpha Ct Nwt Trad D powerCurve powerCurveDomain;
    global u0 windDirections windDistribution Nwd Nws;
    
    %Variable representing whether we are performing the cost with respect
    %to a cost modifier; in this function we want the true cost, i.e.: the
    %expected anual kilowatt per hour.
    applyCostModification = 0;
    
    %The onset wind of each turbine along each wind direction for ever wind
    %slice. For each wind direction w U(1:Nwt,Nws*(w-1)+1:w*Nws) is a Nwt
    %by Nws matrix representing the onset wind of each turbine for each
    %wind slice
    U = zeros(Nwt,Nwd*Nws);
    
    expectedPower=0;

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
        [Q,~] = calcVelocityDeficit( cwD, dwD, inWake, wakeRadius, Trad, applyCostModification);

        %Note that the onset wind in the wind direction w is indexed with
        %respect to the indexed-Turbines.        
        U(:,Nws*(w-1)+1:w*Nws) = calcOnsetWind( Q, dwD, u0, Nwt, Nws, D, alpha, Ct );
        
        %calculate the cost of the wind farm
        turbinePower = interp1(powerCurveDomain, powerCurve, U(:,(w-1)*Nws+1:w*Nws),'linear');
        expectedPower = expectedPower + sum(sum(bsxfun(@times, windDistribution(w,:),turbinePower)));
    end
    
end

%
% end costEvaluation
%