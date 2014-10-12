%{
   We calculate the velocity deficit proportion of turbine i (Ti) shadowing
   turbine j (Tj).

   The (i,j)th entry of Q is
    = 0 if i=j  Since a turbine does nnot shadow itself
    = 0 if i>j  Since by the indexing of the turbines, and by the
                orientation the grid (of the windfarm), and that wind is 
                always comming from the north, we know that Ti is downwind
                to Tj.  Hence the wake of Ti can not act upon Tj.
    = 1 if Tj is in the wake of Ti, and cwD_ij <= wakeRadius_ij - Trad
        i.e.: if Tj is in the wake of Ti and Tj is fully shadowed by Ti
    = X if Tj is in te wake of Ti and Tj is only partially shadowed by Ti's

    We read the value at Q(i,j) as Ti shadows Tj with velocity
    proportion Q(i,j)

   Notes: In the paper "Wind farm micro-sitting by Gaussian particle swarm
   optimization with local search strategy", the authors define the
   velocity deficit of Ti acting on Tj (symbol: q_ij) as:

     q_ij = (A_ij)/(A0)  where A0 = pi*Trad^2
     and
         A_ij = pi*Trad^2 if Tj is fully shadowed by Ti
              =  * if Tj is partially shadowed by Ti
   Now, why not simplify things, since q_ij=1 if Tj is fully shadowed by
   Ti; this is what I did.  Indeed, because of numerical rounding issues,
   (pi*Trad^2)/(pi*Trad^2) does not equal 1.  Hence there will be some
   numerical loss.
   Note: see figure 2 of this paper for an image of full and partial
   shadowing.
  

   Inputs: 
    cwD: (i,j)th entry represents the crosswind distance from Ti to Tj
    dwD: (i,j)th entry represents the downwind distance from Ti to Tj
    wakeRadius: (i,j)th entry represents the radius of the wake from Ti at
                 the location of Tj
    inWake: (i,j)th entry represents whether Tj is in the wake produced
                 by Ti
    Trad: constant turbine-parameter; the radius of the (Vestras v80)
           turbine.

   Returns:
     Q: a Nwt by Nwt matrix whose (i,j)th entry is the velocity deficit of
        Ti shadowing (i.e.: acting on) Tj

Implementation Notes:
   * At the moment, all matrices inputted into this function are of
     size Nwt by Nwt, where Nwt is the number of turbines on the farm.
   * We are using liner indexing with the 'find' command.  By doing so, it
     is quicker to determine whether any turbines are fully or partially
     shadowed. For example, suppose that no turbine is fully shadowed
     (i.e.: no turbine lies fully in another turbine's wake)
      Using matricies:
         shadowedFullySparsity = (cwD<=(wakeRadius-turbRad) & inWake);
      Now, sine no turbine fully lies in another turbine's wake, we have
      shadowedFullySparsity is a (Nwt by Nwt) matrix of all zeros, so to
      make sure we are not indexing by a zero matrix, we must have the
      following code:
         if(any(shadowedFullySparsity(:))~=0) %check for zero matrix
            Q(shadowedFullySparsity)=1;
         end
      Hence, we had to run though the entire matrix and see that every
      element of shadowedFullySparsity is zero.
      Instead, consider using the 'find' command:
         shadowedFullySparsity = find(cwD<=(wakeRadius-turbRad) & inWake);
      Because we are using 'find', shadowedFullySparsity is now an empty
      matrix (not zeros, but empty).  Hence the check to see if
      shadowedFullySparsity is empty is very quick:
         if(any(shadowedFullySparsity)~=0)
            Q(shadowedFullySparsity)=1;
         end
      By using 'find' we avoid having to check if every element of the
      matrix 'shadowedFullySparsity' is zero.
      
%}
function [ Q, shadowedFullySparsity ] = calcVelocityDeficit ...
                           ( cwD, dwD, inWake, wakeRadius, Trad, mode)
                       
% if mode == 1, modify the cost function to move away from the saddle point
% due to the full shadowing.
                       
    Q = zeros(size(dwD));
    
    %Find every (i,j) location where Tj is fully shadowed by Ti
    shadowedFullySparsity = find(cwD<=(wakeRadius-Trad) & inWake);
    
    if( isempty(shadowedFullySparsity)~=1 ) %check for empty matrix
        % add panety to cwD to avoid flat gradient at fullly shadowed positions. 
        % when cwD = wakeRadius-Trad, Q = 1; 
        % when cwD = 0; Q = 2;
        % Q increases as cwD decreases 
        if mode == 1;
% quadratically decreases Q w.r.t. Cwd
%             slope = 5;
%             Q(shadowedFullySparsity) = ... 
%                 1 + slope*(1 - cwD(shadowedFullySparsity)./...
%                 (wakeRadius(shadowedFullySparsity)-Trad) ).^2;
% % linearlly decreases Q w.r.t. Cwd
            slope = 4;
            Q(shadowedFullySparsity) = ... 
                1 - slope*(cwD(shadowedFullySparsity)./...
                (wakeRadius(shadowedFullySparsity)-Trad) - 1);
% % % fixed slope May not be a good idea
%             slope = -30;
%             Q(shadowedFullySparsity) = slope*( cwD(shadowedFullySparsity) -...
%                 wakeRadius(shadowedFullySparsity) + Trad) + 1;
% works for N=30 case
%             C_flat = 1;% constant that controls the slope of the gradient
%             Q(shadowedFullySparsity)=...
%                 (1+C_flat)*(wakeRadius(shadowedFullySparsity)-Trad)./...
%                 (cwD(shadowedFullySparsity)+C_flat*(wakeRadius(shadowedFullySparsity)-Trad));
        else
            Q(shadowedFullySparsity)=1;
        end
    end

    %Find every (i,j) location where Tj is partially shadowed by Ti
    shadowedPartiallySparsity = ...
        find(( (wakeRadius - Trad)<=cwD) & ...
               (cwD<=(wakeRadius + Trad) & inWake));
           
    %if no turbine is partially shadowed, we have finished calculating Q
    %and may leave the function because Q was initialized to zeros

    if( isempty(shadowedPartiallySparsity) ) 
        return;
    end

    %PostCondition: there exists turbines that are partially shadowed, and
    %the variable shadowedPartiallySparsity is a vector of linear indicies
    %of the partially shadowed turbines

    I = shadowedPartiallySparsity;
    
    A0 = pi*Trad^2; %turbine swept area
    alpha1 = 2*acos((Trad^2 + cwD(I).^2 - (wakeRadius(I).^2))./(2*Trad*cwD(I)));
    alpha2 = 2*acos((wakeRadius(I).^2 + cwD(I).^2 - (Trad^2))./(2*wakeRadius(I).*cwD(I)));
    
    Q(I) = (1/(2*A0))*(Trad^2*(alpha1-sin(alpha1)) + (wakeRadius(I).^2).*(alpha2-sin(alpha2)));
%    Q = triu(Q,1); %?????????
end