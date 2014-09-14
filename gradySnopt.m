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
  index4relD: sparsity pattern for the Jacobian G of the constraints F
  u0: freedom wind speed; currently only 12 m/s
  alpha: wake spreading constant
  Ct: turbine thrust constant
  Nwt: number of turbines (constant)
  Trad: radius of wind turbine (constant)
  D: downwind rotor diameter (constant)
%}
global index4relD index4G scaleRelD scaleConstraints4U;
global u0 alpha Ct gamma Nwt Trad D CostModification;
global movements;
%{
    windDirections is an array which specifies which wind directions on a
    rose map we are considering.  Each element in the array must be given
    in degreesm and follow the orientation given in a Rose Map (eg: 0
    degrees is due-north, whereas it would be 90 degrees on a standard
    cartesian coordinate system)
%}
global windDirections Nwd;

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
gamma = 1-sqrt(1-Ct); %this is for the calculation of tildeU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the [SNOPT constraints] with respect to GRADY CASE I
scaleRelD = 1e-6;
scaleConstraints4U = 10;  %this is new; used to be 10
L = 1800;         %windfarm length of side to square
dMin = 4*Tdiam;   %minimum safe (relative) distance between turbines
MinSquaredRelativeDistance = dMin^2;
MaxSquaredRelativeDistance = 2*L^2;
%%%%%%%%%%%%%%%%%%%%%

usrfun = 'gradyUserFun';

u0=12; %units: m/s
windDirections = 0;
Nwd = numel(windDirections);
Nwt = 30;

CostModification = 1;% use modified cost function in SNOPT

%{
   First, we generate an initial guess for Snopt and then use the expected
   power of that initial guess as a lower bound for the cost function (see
   Flow Fupp below.
   Format of xInit:
     -- xInit(1:Nwt)  xPositions of turbines 1 to Nwt
     -- xInit(Nwt+1:2*Nwt) yPositions of turbines 1 to Nwt
     -- x((1+w)*Nwt+1):(2+w)*Nwt) is the onset wind in the w-th wind diretion
          where w=1,...,Nwd
%}
[xInit, expectedPowerInitial ] = buildInitialGuess( Nwt, Nwd, L, dMin, Trad, alpha );

%for making a movie on how Snopt converges to the optimal opsition
movements(:,1) = xInit(1:2*Nwt);

%As there are Nwt turbines and Nwd wind directions: there are Nwt xPositions,
%Nwt yPositions, and Nwt*Nwd Onset wind variables
numDecVar = length(xInit);

%Since there are Nwt turbines on the farm, there are Nwt xPositions and Nwt
%yPositions
numXpos = Nwt;
numYpos = Nwt;
%For each wind direction (Nwd in total) we have Nwt onset wind constraints
%(one for each turbine).  Hence the total number of onset wind desision
%variables we have are Nwd*Nwt.
numOnsetWindVar = Nwd*Nwt;

 %{
   Next we set the lower and upper bounds on the decision variables 
   *Since the positions of the turbines are in [0,L]x[0,L], we have that:
     - Lower bound >=0
     - Upper bound <= L
   *The lower and upper bounds for the onset wind of each turbine are 0  and
    u0 respectivly 
 %}
xLow = zeros(numDecVar,1); 
xUpp = [ L*ones(numXpos,1);    %upper bound for xPositions
         L*ones(numYpos,1);    %upper bound for yPositions
         u0*ones(numOnsetWindVar,1) ]; %upper bound for the onset wind of each turbine for each wind direction


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



% SparsityPattern: to be finished
jacobianSparsityConstraints = findSparsityPattern(Nwt,Nwd,index4relD);
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
perfectPower = (0.3*(u0^3*Nwt))*Nwd;
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