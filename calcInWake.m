%{
   This function calculates and return a matrix  whose the (i,j)th entry is
             = 1 if turbine j is in the wake of turbine i
             = 0 otherwise.

   Turbine j (Tj) is in the wake produced by turbine i  (Ti) if the
   downwind distance between Ti and Tj is greater than zero, and if the
   crosswind distance between Ti and Tj is less than the radius of the wake
   from Ti at the location of Tj plus the radius of the turbine itself.

   In math, we have that Tj is in the wake of Ti if dwD(i,j)>0 and
   cwD(i,j) < wakeRadius(i,j) +Trad

   Inputs:
     cwD  matrix of cross wind distances between turbines
   
   Output:
     inWake  a Logical matrix whose the (i,j)th entry is
                  = 1 if turbine j is in the wake of turbine i
                  = 0 otherwise.

   Implementation Notes:
      We do not check if the downwind distance between Ti and Tj is greater
      than zero.  If dwD(i,j)=0 then we know that cwD(i,j)<wakeRadius(i,j)
      +Trad  will fail.  It will fail because dwD(i,j)=0 then we know that
      wakeRadius(i,j)=turbineRadius, and because of the Minimum distance
      safty requirements of the turbines, we will never this situation
   Notes: 
     *inWake has the same sparsity pattern as dwD
%}
function [ inWake ] = calcInWake( cwD, wakeRadius, TurbRad )
   inWake = cwD < wakeRadius + TurbRad*triu(ones(size(cwD)),1);
end

