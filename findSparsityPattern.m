%{
         | x| y |U_1|...|U_Nwd|
      ---------------------------
   cost  |  |   |   |   |     |
   relD  |  |   |   |   |     |
   U_1   |  |   |   |   |     |
   ...  
   U_Nwd |  |   |   |   |     |
%}
function [ jacobianSparsityConstraints ] = findSparsityPattern( Nwt, Nwd, index4relD )

%first turn the logical indexing into equivalent row and column indicies
[sparsityRow, sparsityCol] = ind2sub([Nwt,Nwt],index4relD);

%{
   Here we build the sparsity pattern for the relative distance
   constraints for the jacobian of F (which contains the constraints to the
   problem)
   Notes:
     * numel(sparsity)=numel(sparsityRow)=numel(sparsityCol)
%}
numRelDcon = numel(index4relD);
relDjacobiSparsity = zeros(numRelDcon,Nwt);
for i=1:numRelDcon
    relDjacobiSparsity(i,sparsityRow(i))=1;
    relDjacobiSparsity(i,sparsityCol(i))=1;
end

%{
   jacobianSparsity represents the (sparsity) of the jacobian of F defined
   in gradySnopt.m; in gradySnopt.m F is built with the cost function on
   row one, and then the (squared) relative distance constraints.  We make
   the cost function dense: that is, we make SNOPT calculate all jacobian
   entries for the cost function

   Note: the decision variables are [x(1),..., x(Nwt), y(1),..., y(Nwt)]
         where [x(i),y(i)] is the position of the ith turbine
%}

%%%
%%% The following are multiple versions of jacobian sparsity
%%%

%---sparsity pattern for cost function and relative distances---%
%{
   This just builds a block diagonal matrix of lower triangular matricies
   along the main diagonal (of matricies)
   Note: the '+2' is from Nwd columns, along with one col for x and one col
   for y, each with Nwt entries
%}
sparsOnsetWindCnstrWrtU = zeros(Nwd*Nwt,Nwd*Nwt);
for w=1:Nwd
    sparsOnsetWindCnstrWrtU((w-1)*Nwt+1:w*Nwt, (w-1)*Nwt+1:w*Nwt) = tril(ones(Nwt,Nwt));
end

jacobianSparsityConstraints = [ zeros(1,2*Nwt), ones(1,Nwd*Nwt) ; ...
                    horzcat(relDjacobiSparsity,relDjacobiSparsity,zeros(numRelDcon,Nwd*Nwt));
                    horzcat(ones(Nwd*Nwt,2*Nwt),sparsOnsetWindCnstrWrtU) ];

end

%{

%---dese with respect to the onset wind constraints
% jacobianSparsityConstraints = [ zeros(1,2*Nwt), ones(1,Nwd*Nwt) ; ...
%                     horzcat(relDjacobiSparsity,relDjacobiSparsity,zeros(numRelDcon,Nwd*Nwt));
%                     ones(Nwd*Nwt,Nwt*(Nwd+2))];
                

%---no spars patt for cost fun---%
% jacobianSparsityConstraints = [ ones(1,2*Nwt), ones(1,Nwt) ; ...
%                     horzcat(relDjacobiSparsity,relDjacobiSparsity,zeros(numRelDcon,Nwt));
%                     ones(Nwt,3*Nwt)];

%---Fully dense: all ones
% jacobianSparsityConstraints = [ ones(1,2*Nwt), ones(1,Nwt) ; ...
%                     ones(numRelDcon,3*Nwt);
%                     ones(Nwt,3*Nwt)];
                 
%}

%
% end findSparsityPattern.m
%