%The initial guess for the turbine positions is the three-row turbine
%positions

%the turbine positions (which consists of one x and one y value), and onset wind speed are our
%decision variables

% Note: the decision variables are [x(1),..., x(Nwt), y(1),..., y(Nwt), U(1), ... U(Nwt)]
%          where [x(i),y(i)] is the position of the ith turbine and U(i)
%          are the onset wind speed


clear all;
%close all;

%{
  Index_4_relative_D: sparsity pattern for the Jacobian G of the constraints F
  u0: freedom wind speed; currently only 12 m/s
  alpha: wake spreading constant
  Ct: turbine thrust constant
  Nwt: number of turbines (constant)
  Trad: radius of wind turbine (constant)
  D: downwind rotor diameter (constant)
%}
global Index_4_relative_D Index_4_G scale_Rel_D;
global u0 alpha Ct gamma Nwt Trad D CostModification;

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
scale_Rel_D = 1e-6;
L = 1800;         %windfarm length of side to square
dMin = 4*Tdiam;   %minimum safe (relative) distance between turbines
MinSquaredRelativeDistance = dMin^2;
MaxSquaredRelativeDistance = 2*L^2;
%%%%%%%%%%%%%%%%%%%%%

u0=12; %units: m/s

usrfun = 'gradyUserFun';

%Make an initial guess for SNOPT
%{
  Here we can use the functions 'gradyScheme' or 'empiricalScheme' to make
  a nice guess.  Everything is in column vectors.
%}

Nwt = 5;

%%%Initial Guess
turbineSites = 1800*rand(2*Nwt,1); %---- random turbine placement
[expectedPower, U] = Cost_evaluation( turbineSites );
turbineSites = [ turbineSites(1:Nwt), turbineSites(Nwt+1:2*Nwt), U ];
XLOC = 1;  %enumeration: location of xPositions in turbineSites is col 1
YLOC = 2;  %enumeration: location of yPositions in turbineSites is col 2
ULOC = 3;  %enumeration: location of onset winds for turbines is col 3
indexedTurbines = sortrows(turbineSites,[-YLOC,XLOC]);
turbineXpos = indexedTurbines(:,XLOC); %loc of x values of turbines in col vec
turbineYpos = indexedTurbines(:,YLOC); %loc of y vals of turbs in col vec
Onsetwind   = indexedTurbines(:,ULOC); % onset wind speed 
xInit = [turbineXpos; turbineYpos; Onsetwind]; %xInit is a column vector
%%%end initial guess

numDecVar = length(xInit);

%set lower and upper bounnds on the decision variables
 %{
    Since the positions of the turbines are in [0,L]x[0,L], we have that:
      * Lower bound >=0
      * Upper bound <= L
 %}
xLow = zeros(numDecVar,1); 
xUpp = [L*ones(2*Nwt,1); u0*ones(Nwt,1)];


%determine the number of relative-distance constraints for the problem.
%For every two turbines there must exist a minimum (saftey) distance
%between them.  'sparsity' represents the locations where we need to
%enforce the relative distance constraint: it's an upper triangular matrix
%with zeros along the main diagonal. since relative distance is symmetric.
%i.e.: the distance between T1 and T2 is the same as the distance between
%T2 and T1
Index_4_relative_D = find(triu(ones(Nwt,Nwt),1));
numRelDistConstraints = numel(Index_4_relative_D);

%{
   Upper and Lower Constraints:
     1) Cost Function (the expected power of a wind farm configuration)
     2) Relative (safe) minimum distance between each turbine
     3) Constraints that define onset wind speed U
%}
Flow = [ expectedPower;
         scale_Rel_D*MinSquaredRelativeDistance*ones(numRelDistConstraints,1);
         zeros(Nwt,1)]; 
         
Fupp = [ 1.0e+10;
         scale_Rel_D*MaxSquaredRelativeDistance*ones(numRelDistConstraints,1);
         zeros(Nwt,1)];
     
%number of constraints including the cost function   
numConstraints = length(Flow); 

xMul   = zeros(numDecVar,1); xState = zeros(numDecVar,1);
Fmul   = zeros(numConstraints,1); Fstate = zeros(numConstraints,1);
ObjAdd = 0; ObjRow = 1;
A  = [];  iAfun = [];  jAvar = [];

% SparsityPattern: to be finished
jacobianSparsityConstraints = findSparsityPattern(Nwt,Index_4_relative_D);
[iGfun,jGvar] = find(jacobianSparsityConstraints); 
Index_4_G = find(jacobianSparsityConstraints); 

%Optimal Parameters for SNOPT. See chapter 7, pg 62 'Optimal Parameters'
%Note we first set 'Defaults' to start SNOPT on a clean-slate; very important!
snset ('Defaults'); %<-- Indeed, without this SNOPT could give false results
snseti( 'Major Iteration limit', 5000); % SNOPT userg guide pp. 62 
snseti('Derivative option', 0);
snseti('Verify level', 3);
snset('Maximize');
primal_tol =1*1.0e-6;
snsetr   ( 'Major feasibility tolerance', primal_tol);
snsetr   ( 'Minor feasibility tolerance', primal_tol);
dual_tol = 1*1.0e-6;
snsetr   ( 'Major optimality tolerance',  dual_tol);
snsetr   ( 'Minor optimality tolerance',  dual_tol);
snsetr   ( 'Major print level', 3);


solveopt = 1;
numberRuns = 100;
perfectPower = 0.3*u0^3*Nwt;
tolPercent = 1; %if 

xOptsWithCostMod = zeros(2*Nwt, numberRuns);
xOptsNoCostMod = zeros(2*Nwt, numberRuns);

expPwrUsingCostMod = zeros(numberRuns, 1);
expPwrNoCostMod = zeros(numberRuns, 1);

for i=1:numberRuns

    CostModification = 1;
    [xOptWithCostMod,~,~,~,~] = snoptcmex( solveopt, ...
                                         xInit,xLow,xUpp,xMul,xState, ...
                                         Flow,Fupp,Fmul,Fstate,        ...
                                         ObjAdd,ObjRow,A,iAfun,jAvar,  ...
                                         iGfun,jGvar,usrfun );

    powerUsingCostMod = Cost_evaluation( xOptWithCostMod(1:2*Nwt) );
    
    CostModification = 0;
    [xOptNoCostMod,~,~,~,~] = snoptcmex( solveopt, ...
                                         xInit,xLow,xUpp,xMul,xState, ...
                                         Flow,Fupp,Fmul,Fstate,        ...
                                         ObjAdd,ObjRow,A,iAfun,jAvar,  ...
                                         iGfun,jGvar,usrfun );

    powerWithoutCostMod = Cost_evaluation( xOptNoCostMod(1:2*Nwt) );

    
    xOptsWithCostMod(:,i)=xOptWithCostMod(1:2*Nwt) ;
    xOptsNoCostMod(:,i)=xOptNoCostMod(1:2*Nwt) ;
    expPwrUsingCostMod(i)=powerUsingCostMod;
    expPwrNoCostMod(i)=powerWithoutCostMod;
    
    if(mod(i,50)==0)
        disp(num2str(i));
    end
    
    
    %Make new Initial Guess
    turbineSites = L*rand(2*Nwt,1);
    [expectedPower, U] = Cost_evaluation( turbineSites );
    turbineSites = [ turbineSites(1:Nwt), turbineSites(Nwt+1:2*Nwt), U ];
    XLOC = 1;  %enumeration: location of xPositions in turbineSites is col 1
    YLOC = 2;  %enumeration: location of yPositions in turbineSites is col 2
    indexedTurbines = sortrows(turbineSites,[-YLOC,XLOC]);
    turbineXpos = indexedTurbines(:,1); %loc of x values of turbines in col vec
    turbineYpos = indexedTurbines(:,2); %loc of y vals of turbs in col vec
    Onsetwind   = indexedTurbines(:,3); % onset wind speed 
    xInit = [turbineXpos; turbineYpos; Onsetwind]; %xInit is a column vector
    %----end initial guess
    
    
end%-end test run

%Analize results using relative (percentage) error
%Note: a sucsess is when the expected power of the wind farm is within one
%percent of the perfecet power
amountWithinPerfectUsingCostMod =  ...
         numel(find( abs((expPwrUsingCostMod-perfectPower)/perfectPower)*100<tolPercent));
     
amountWithinPerfectNoCostMod = ...
         numel(find( abs((expPwrNoCostMod-perfectPower)/perfectPower)*100<tolPercent));

percentageSuccessWithCostMod = amountWithinPerfectUsingCostMod/numberRuns*100;     
percentageSuccessNoCostMod = amountWithinPerfectNoCostMod/numberRuns*100;
          
fileName = strcat('costModComparison_Nwt_',num2str(Nwt));

save(fileName, 'xOptsWithCostMod', ...
               'xOptsNoCostMod', ...
               'expPwrUsingCostMod', ...
               'expPwrNoCostMod', ...
               'amountWithinPerfectUsingCostMod',...
               'amountWithinPerfectNoCostMod', ...
               'percentageSuccessWithCostMod',...
               'percentageSuccessNoCostMod', ...
               'Nwt', 'numberRuns' );
