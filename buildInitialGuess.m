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
function [ xInit, expectedPower ] = buildInitialGuess( Nwt, Nwd, Nws, L, ...
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


xInit = [turbineSites(1:Nwt); turbineSites(Nwt+1:2*Nwt); U(:)];

end

%
% end buildInitialGuess.m
%