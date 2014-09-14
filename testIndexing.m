%{
 Question:
   When Snopt moves the turbines about to search for an optimal point, does
   it break the left-right top down indexing (which is needed to calculate
   dwD and cwD)

   We want to know how Snopt moves the turbines about when searching for
   the optimal solution.  At each call to userFun, we need to reindex

 Format:   
   movements: each column is a turbineSite ordered in the standard fashion
   where indicies 1 to Nwt are the xPositions and Nwt+1 to 2*Nwt 

   movements(:,1) == the inital guess (stored in gradySnopt in the initial
                     guess section
   movements(:,end) == the (suposed) optimal solution that Snopt returns,
                       stored at each iteration of the userFun (which is
                       called by Snopt)
%}

Nwt = size(movements,1)/2;
totalMovements = size(movements,2);

%check each turbine movement for a violation of calculation of downwind
%distance.  A voilation occurs if an upwind turbine has been moved downwind
%by Snopt
for i=1:totalMovements
    turbineSites = movements(:,i);  %column vec
    dwD = calcDownWindDist(turbineSites(Nwt+1:2*Nwt)); %dwD needs just y-values
    
    %look for atleast one situation where an element of dwD is less than
    %zero
    foundViolation = any(find(dwD<0));
    if( foundViolation )
        disp(strcat('Violation of dwD at index ', num2str(i) ));
    end
end

if( ~foundViolation )
    disp('Down Wind distance was correctly computed');
end

