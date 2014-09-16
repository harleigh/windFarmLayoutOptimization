%{
   Build an initial guess for Snopt

   Inputs: Various parameters to computer various initial turbine
   placements on the windFarm.
         Nwt: Number of wind turbines
         L:   Length of square wind farm
         dMin: minimum safe distance
         Trad: radius of a turbine
         alpha: wake spreading constant

   Outputs:
     xInit: Initial guess for Snopt, in colum format
            xInit(1:Nwt) == xPositions of turbines
            xInit(Nwt+1:2*Nwt) == yPositions of turbines
            xInit(2*Nwt+1: (Nwd+2)*Nwt) == onset wind of each turbine along
               each wind direction w.  Eg: for the onset wind along
               direction w we have x((w+1)*Nwt: (w+2)*Nwt)

         Example: if Nwt = 30, then [xInit(1), xInit(31), xInit(61)] is the x,y
                  and onset wind of the first turbine

     expectedPower: Expected anual power generation of wind farm in
                    Kilowatts per hour

   Notes:
     The building of an initial guess was removed from the main snopt file
     as by isolating it here, it is easier to tweek the initial guess,
     rather than scroll through the main snopt file.
%}
function [ xInit, expectedPower ] = buildInitialGuess( Nwt, Nwd, L, ...
                                                       dMin, Trad, alpha )

%Generate initial guess for turbine x,y positions on a square farm of width
%L. Note: the x,y position are such that the first Nwt entries are the
%xPositions x1 x2,...,xNwt and the next Nwt entries are the yPositions y1
%y2,..,yNwt, where turbine 1 is at (x1,y1)
    % turbineSites = gradyScheme(1,'');
    % turbineSites = empiricalScheme(L,Nwt);
    % turbineSites = globalScheme(L, Nwt, dMin, Trad,alpha);
%      load paperResult gPop;
%      turbineSites = gPop;
% load 'currentRunResults' 'xOpt'
% turbineSites = xOpt(1:2*Nwt);
  %turbineSites = L*rand(2*Nwt,1);
turbineSites = randomScheme( Nwt, L, dMin );


%next we computer the true cost (no cost modification) of the turbine sites
%that we generated above.  That is determine the expected power produced
%for this given configuration of turbines and determine the onset wind to
%each turbine for each wind direction.
[expectedPower, U] = trueCostEvaluation( turbineSites );

turbineProfiles = [ turbineSites(1:Nwt), turbineSites(Nwt+1:2*Nwt), U ];
XLOC = 1;  %enumeration: location of xPositions in turbineProfiles is col 1
YLOC = 2;  %enumeration: location of yPositions in turbineProfiles is col 2
%ULOC is the location of the onset wind constraints for each wind
%direction begins at column 3 and has Nwt columns (each of length Nwt); the '+2'
%is because the first two columns are x,y position
ULOC = 3:(Nwd+2);
turbineXpos = turbineProfiles(:,XLOC); %loc of x values of turbines in col vec
turbineYpos = turbineProfiles(:,YLOC); %loc of y vals of turbs in col vec
Onsetwind   = turbineProfiles(:,ULOC); % onset wind speed of turbines: Nwt by Nwd

%pack the generated guess into Snopt format--i.e.: as a clolumn
xInit = [turbineXpos; turbineYpos; Onsetwind(:)];

end

%
% end buildInitialGuess.m
%