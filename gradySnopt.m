%{
This version of optimizing wind farm turbine placement.
Onset wind is a decision variable
Optimizes with respect multiple wind directions at one speed
Contains multiple cost modifications (see calcVelocityDeficit)
Contains multiple sparsity Patterns for the Jacobian (see
  findSparsityPattern)
Scales constraints to prioitize Snopts decision on which constraint is most
  important
%}

%The initial guess for the turbine positions is the three-row turbine
%positions

%the turbine positions (which consists of one x and one y value), and onset wind speed are our
%decision variables

% Note: the decision variables are [x(1),..., x(Nwt), y(1),..., y(Nwt), U(1), ... U(Nwt)]
%          where [x(i),y(i)] is the position of the ith turbine and U(i)
%          are the onset wind speed


clear all;
close all;

%{
    Due to SNOPT, we must use global variables for any data that we wish to
    use in the userfun.  Globals are used only in two places:
      *in the userfun
      *trueCostEvaluation (just to avoid a massive amount of parameters)

  ---SNOPT (special) parameters---
  index4relD: sparsity pattern for the Jacobian G of the constraints F
  index4G: a logical index identifying the nonzero elements of G
  scaleRelD: scale-factor to (squared) relative distance between turbines
  scaleConstraints4U: scale factor for the constraints on the onset wind U
  CostModification: determines whether we are running a modifed cost
                    function (used to help avoid local minima)

  ---Turbine Parameters---
  alpha: wake spreading constant
  Ct: turbine thrust constant
  Nwt: number of turbines (constant)
  Trad: radius of wind turbine (constant)
  D: downwind rotor diameter (constant)
  powerCurve: the power curve of the turbine used in the simulation
              (currently Vestras-80)
  powerCurveDomain: wind speeds that the power curve are defined upon; used
                    for interpolation purposes

  movements: this exists to capture how SNOPT converges to the 'optimal'
             solution (where later we make a movie with it)

  ----Wind Parameters---
  u0: a vector, free wind speed
  Nws: Number of wind slices (amount that u0 is discritized; num of elements in u0)
  windDirections: the directions of wind active in the rose map
  Nwd: Number of wind directionss
  windDistribution: with windDirections, this constitutes the data of a
                    rose map
%}
global index4relD index4G scaleRelD scaleConstraints4U CostModification;
global alpha Ct gamma Nwt Trad D powerCurve powerCurveDomain;
global movements;
global u0 windDirections windDistribution Nwd Nws;


%First we set the parameters of the simulation: turbine wind and Snopt
%paramaters.
[alpha, Ct, gamma, Nwt, Trad, D, L, dMin, ...
     powerCurve, powerCurveDomain] = setTurbineParameters();
 
[u0, windDirections, windDistribution, Nwd, Nws] = setWindParameters();
%now set the snopt parameters 
CostModification = 1;%whether we want to apply the modified cost calculation (to avoid local optima)
scaleRelD = 1e-6;
scaleConstraints4U = 10;
MinSquaredRelativeDistance = dMin^2;
MaxSquaredRelativeDistance = 2*L^2;
%%%%%%%%%%%%%%%%%%%%%----end setting parameters




%userfun contains the calculation of the exected power of a turbine
%setting, it packs the constraints into F, and calculates as much as G (the
%jacobian of F wrt the decision variables) as possible.
usrfun = 'gradyUserFun';

%{
   First, we generate an initial guess for Snopt and then use the expected
   power of that initial guess as a lower bound for the cost function (see
   Flow Fupp below).
   Format of xInit:
     -- xInit(1:Nwt)  xPositions of turbines
     -- xInit(Nwt+1:2*Nwt) yPositions of turbines
     -- x((2+(w-1)*Nws)*Nwt+1:(2+w*Nws)*Nwt) is the onset wind in the w-th wind diretion
          where w=1,...,Nwd
%}
[xInit, expectedPowerInitial ] = buildInitialGuess( Nwt, Nwd, Nws, ... 
                                                    L, dMin, Trad, alpha );

%for making a movie on how Snopt converges to the optimal opsition
movements(:,1) = xInit(1:2*Nwt);

%This is used to convienently set the upper and lower bounds of the
%decision variables (see xLow and xUpp)
numDecVar = length(xInit);

%There are Nwt turbines on the 2D farm, hence there are Nwt xPositions as well as
%Nwt yPositions
numXpos = Nwt;
numYpos = Nwt;
%For each wind direction (Nwd in total) we have Nwt*Nws onset wind constraints
%(one for each turbine, for each wind slice).  Hence the total number of onset
%wind desision variables we have are Nwd*Nwt*Nws.
numOnsetWindVar = Nwd*Nwt*Nws;

 %{
   Next we set the lower and upper bounds on the decision variables.
   *Since the positions of the turbines are in [0,L]x[0,L], we have that:
     - Lower bound of turbine (x,y) positions >=0
     - Upper bound of turbine (x,y) positions <= L
   *The lower and upper bounds for the onset wind of each turbine are 0  and
    the max element in u0 respectivly 
 %}
xLow = zeros(numDecVar,1); 
xUpp = [ L*ones(numXpos,1);    %upper bound for xPositions
         L*ones(numYpos,1);    %upper bound for yPositions
         max(u0)*ones(numOnsetWindVar,1) ]; %upper bound for the onset wind of each turbine for each wind direction


%determine the number of relative-distance constraints for the problem.
%For every two turbines there must exist a minimum (saftey) distance
%between them.  'index4relD' represents the locations where we need to
%enforce the relative distance constraint: it's an upper triangular matrix
%with zeros along the main diagonal. since relative distance is symmetric.
%i.e.: the distance between T1 and T2 is the same as the distance between
%T2 and T1
index4relD = find(triu(ones(Nwt,Nwt),1));
numRelDistConstraints = numel(index4relD);

%each onset wind constraint uses the onset wind, so the number of onset
%wind constraints is the same as the number of onset wind desision
%variables
numOnsetWindConstraints = numOnsetWindVar;

%{
   Upper and Lower bounds of the Constraints:
     1) Cost Function, the expected power of a wind farm configuration has
        a lower bound of the expected power from the initial guess of the
        turbine replacement
     2) Relative (safe) minimum distance between each turbine.  Note that
        the values of this constraint are really large when compared to the
        constraints of the onset wind, hence we must introduce a scaling
        factor to bring the values down.
     3) The Constraints that define onset wind speed U is an equality
        constraint hence zero
%}
Flow = [ expectedPowerInitial;
         scaleRelD*MinSquaredRelativeDistance*ones(numRelDistConstraints,1);
         zeros(numOnsetWindConstraints,1)]; 
Fupp = [ 1.0e+10;
         scaleRelD*MaxSquaredRelativeDistance*ones(numRelDistConstraints,1);
         zeros(numOnsetWindConstraints,1)];
numConstraints = length(Flow); 



xMul   = zeros(numDecVar,1); xState = zeros(numDecVar,1);
Fmul   = zeros(numConstraints,1); Fstate = zeros(numConstraints,1);
ObjAdd = 0; ObjRow = 1;
A  = [];  iAfun = [];  jAvar = [];



%Next we determine the sparsity of the Jacobian of the constraints with
%respect to the decision variables.
jacobianSparsityConstraints = findSparsityPattern(Nwt,Nwd, Nws, index4relD);
[iGfun,jGvar] = find(jacobianSparsityConstraints); 
index4G = find(jacobianSparsityConstraints); 


%Optimal Parameters for SNOPT. See chapter 7, pg 62 'Optimal Parameters'
%Note we first set 'Defaults' to start SNOPT on a clean-slate; very important!
snset ('Defaults'); %<-- Indeed, without this SNOPT could give false results
snseti('Major Iteration limit', 5000); % SNOPT userg guide pp. 62 
snseti('Derivative option', 0);
snseti('Verify level', 3);
snset('Maximize');
primal_tol =1*1.0e-6;
    snsetr('Major feasibility tolerance', primal_tol);
    snsetr('Minor feasibility tolerance', primal_tol);
dual_tol = 1*1.0e-6;
    snsetr('Major optimality tolerance',  dual_tol);
    snsetr('Minor optimality tolerance',  dual_tol);
snsetr('Major print level', 3);


%Sumary and Solution files; see chapter 8 of SNOPT guide (section 8.8, 8.9)
snprint('resultGradySNOPT.txt'); %print detained info on the progress of Snopt to file
snsummary('resultGradySNOPT.sum'); %print summerized info on the progess of Snopt to file


%call snopt
solveopt = 1;
tic
[xOpt,F,xMul,Fmul,INFO] = snoptcmex( solveopt, ...
 				                     xInit,xLow,xUpp,xMul,xState, ...
 				                     Flow,Fupp,Fmul,Fstate,        ...
 				                     ObjAdd,ObjRow,A,iAfun,jAvar,  ...
 				                     iGfun,jGvar,usrfun );
runTime=toc;
%%%%%%%%%%%%%%%%%%%%%% Display and Plot Code %%%%%%%%%%%%%%%%%%%%%%

snsummary off; % close the summary .sum file; empties the print buffer
snprint   off; % Closes the print-file .out; empties the print buffer

%---display SNOPT results to command line---%
if CostModification == 1
    disp(strcat('With cost modification'));
else
    disp(strcat('Without cost modification'));
end

disp(strcat('Number of turbines: Nwt==', num2str(Nwt)));
disp(strcat('Wind Direction(s): ', num2str(windDirections)));
disp(strcat('Free Wind Speed: u0==', num2str(u0)));

disp(strcat('Execution time: ', num2str(runTime)));
disp(strcat('SNOPT exited with INFO==', num2str(INFO)));

%To calcuate the efficiency of the wind-turbine layout.  Note: Perfect
%Power along Nwd wind directions is equal to perfect power along each wind
%direction.  
perfectPower = calcPerfectPower(Nwt, Nwd, u0, windDistribution, powerCurveDomain, powerCurve);
[expectedPowerOpt,U] = trueCostEvaluation( xOpt(1:2*Nwt) );
disp(strcat('Expected Power of Initial Guess: ', num2str(expectedPowerInitial)));
disp(strcat('-----Efficiency: ', strcat(num2str(expectedPowerInitial/perfectPower*100),'%')));
disp(strcat('Expected Power of Optimal: ', num2str(expectedPowerOpt)));
disp(strcat('-----Efficiency: ', strcat(num2str(expectedPowerOpt/perfectPower*100),'%')));
fprintf('\n');
%When plotting the turbine sites, Note that indexing (1:2*Nwt) grabs the x then
%the y positions of the turbines.
showWake=0;
includeTurbineRadius=1;
plotTurbineScheme(xInit(1:2*Nwt), L, showWake, [], includeTurbineRadius, 'Initial Guess' );
plotTurbineScheme(xOpt(1:2*Nwt), L, showWake, [], includeTurbineRadius, 'Optimal Solution');

save currentRunResults;
%
% end gradySnopt.m
%

%{
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parameters for the Vestas V80-1800 wind turbine for GRADY CASE I
% h = 60;                %hub height
% Tdiam = 40;            %diameter of turbine blades
% Trad = Tdiam/2;        %radius of turbine blades
% z0 = 0.3;              %roughness constant
% alpha = 0.5/log(h/z0); %wake spreading constant
% Ct = 0.88;             %turbine thrust coefficient
% a=(1-sqrt(1-Ct))/2;    %for calculation of D on next line
% D = Tdiam*sqrt((1-a)/(1-2*a)); %downstream rotor diameter
% gamma = 1-sqrt(1-Ct); %this is for the calculation of tildeU
% 
% Nwt = 10;
% powerCurve = [0,0,0,0,2,97,255,459,726,1004,1330,1627,1772,1797,1802 * ones(1,10),1800,1800,zeros(1,5)];
% powerCurveDomain = [0:30];
% L = 1800;         %windfarm length of side to square
% dMin = 4*Tdiam;   %minimum safe (relative) distance between turbines

%------Wind Parameters
% u0=[8 12 17]; %units: m/s, u0 must be a row vector
% Nws = numel(u0);
% windDirections = [0 10 20];
% Nwd = numel(windDirections);
% 
% %windDistribution must be Nwd by Nws matrix
% windDistribution = [ones(2,1)*[0.0039    0.0085    0.0115];
%                     0.0039    0.0108    0.0135];
%}