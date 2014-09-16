%{
   Purpose: Given a square wind farm of length L return a random turbine
   placement of Nwt turbines that satisfy the minimum relative safe
   distance (dMin) saftey requirement

   Outputs:
    turbineSites: a column vector of length 2*Nwt containing the (x,y)
    locations of the turbines on the wind farm such that
       -- turbineSites(1:Nwt) are the x positions
       -- turbineSites(Nwt+1:2*Nwt) are the y positions
       -- the relative distance between any two (distinct) turbines is
           atleast dMin
    eg: (turbineSites(i),turbineSites(i+Nwt)) is the (x,y) position of the
        ith turbine i=1,...,Nwt

%}
function [ turbineSites ] = randomScheme( Nwt, L, dMin )

    turbineSites = zeros(2*Nwt,1);
    turbinePos = L*rand(2,1);   %generate (x,y) position in wind farm
    turbineSites(1) = turbinePos(1);        %store x position
    turbineSites(1+Nwt) = turbinePos(2);    %store y position
    
    numberTriesAttempted = 0;
    MaxPlacementAttempts = 10000;
    
    i=2; %iterator for storing turbine postions in turbineSites matrix
    while (i <= Nwt && numberTriesAttempted < MaxPlacementAttempts )
        turbinePos = L*rand(2,1);
        %note: (i-1) is equal to the number of turbines already stored in
        %turbineSites
        validPlacement = isValidPosition(turbinePos, (i-1) ,turbineSites, Nwt, dMin);
        if( validPlacement )
            turbineSites(i) = turbinePos(1);
            turbineSites(i+Nwt) = turbinePos(2);
            i = i+1;
            numberTriesAttempted=0;
        else
            numberTriesAttempted=numberTriesAttempted+1;
        end
    end
    
    %Check if we were not able to place a turbine on the farm; this is
    %rare.
    if ( numberTriesAttempted >= MaxPlacementAttempts )
        disp('Could not satisfy minimum safety requirements on random scheme');
    end
    

end

function [ satisfiesMinimumSafeDist ] =  ...
     isValidPosition( turbinePos, numTurbinesPlaced, turbineSites, Nwt, dMin )
    
    satisfiesMinimumSafeDist = true;
    
    x = turbinePos(1);
    y = turbinePos(2);
    
    for i=1:numTurbinesPlaced
        curX = turbineSites(i);
        curY = turbineSites(i+Nwt);
        relativeDistance = sqrt( (x-curX)^2 + (y-curY)^2 );
        if( relativeDistance < dMin )
            satisfiesMinimumSafeDist = false;
            break;  %break and return 'satisfiesMinimumSafeDist'
        end
    end
    
end

