%{
   Given turbine positions on a square grid of length L, which is oriented
   as:
                   (0,L)    (L,L)
                     *-------*
                     |       |
                     |       |
                     *-------*
                   (0,0)    (0,L)
   Calculate the expected power output of this configuration when the
   wind is 12 m/s ONLY from the north--i.e.: 0 degrees

   Special Notation (in comments):
     * Ti, Tj is turbine i, and turbine j
     * We always assume that all of the turbines on the wind farm are of
       the same type
     * The number of wind turbines on the farm is fixed (i.e.: already
       given)

   Inputs:
     x: A column vector of Nwt turbine positions ordered such that
        x(1:Nwt) are the x-positions of the turbines and x(Nwt+1:2*Nwt) are
        the y-positions of the turbines
  
   Outputs:
     expectedPower: The expected power from the wind farm
%}
function [ expectedPower ] = testOneDirectionOneSpeedCost( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the Vestas V80-1800 wind turbine for GRADY CASE I
h = 60;                %hub height
Tdiam = 40;            %diameter of turbine blades
Trad = Tdiam/2;        %radius of turbine blades
z0 = 0.3;              %roughness constant
alpha = 0.5/log(h/z0); %wake spreading constant
Ct = 0.88;             %turbine thrust coefficient
a=(1-sqrt(1-Ct))/2;    %for calculation of D on next line
D = Tdiam*sqrt((1-a)/(1-2*a)); %downstream rotor diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end Parameters section

%Number of wind turbines.
Nwt = numel(x)/2;

%free wind speed, 12 m/s
u0=12;

turbineSites = [x(1:Nwt) x(Nwt+1:2*Nwt)];

XLOC = 1;  %enumeration: location of xPositions in turbineSites is col 1
YLOC = 2;  %enumeration: location of yPositions in turbineSites is col 2

%TODO: here is where we would rotate for the wind

%index the turbines by downwind distance; by the orientation of the wind
%farm, the first turbines into the wind are the ones with the higest
%y-value
indexedTurbines = sortrows(turbineSites,-YLOC);

%{
  Compute the following constants before the wake calculation:
    dwD: (i,j)th entry represents the downwind distance from Ti to Tj
    cwD: (i,j)th entry represents the crosswind distance from Ti to Tj
    wakeRadius: (i,j)th entry represents the radius of the wake from Ti at
                 the location of Tj
    inWake: (i,j)th entry represents whether Tj is in the wake produced
                 by Ti
    Q: (i,j)th entry represents the velocity deficit of Ti shadowing Tj
%}
dwD = calcDownWindDist(indexedTurbines(:,YLOC));
cwD = calcCrossWindDist(indexedTurbines(:,XLOC));
wakeRadius = calcWakeRadius(dwD, Trad, alpha);
inWake = calcInWake(cwD, wakeRadius);
Q = calcVelocityDeficit( cwD, dwD, inWake, wakeRadius, Trad);
U = calcOnsetWind( Q, dwD, u0, D, alpha, Ct );

%currently the expected power is calculated as below, but we can also
%calculate it using the power curve.
expectedPower = sum(U.^3)*0.3;


end

