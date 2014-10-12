%{
         | x| y |U_1|...|U_Nwd|
      ---------------------------
   cost  |  |   |   |   |     |
   relD  |  |   |   |   |     |
   U_1   |  |   |   |   |     |
   ...  
   U_Nwd |  |   |   |   |     |
%}
function [ jacobianSparsityConstraints ] = findSparsityPattern( Nwt, Nwd, Nws, index4relD )

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

%---sparsity pattern for cost function and relative distances---%
%{
   This just builds a block diagonal matrix of lower triangular matricies
   along the main diagonal (where each entry is a matrix of Nwt by Nwt)
   Note: the '+2' is from Nwd columns, along with one col for x and one col
   for y, each with Nwt entries
%}
sparsityAlongOneWindDir = zeros(Nws*Nwt, Nws*Nwt);
sparsOnsetWindCnstrWrtU = zeros(Nwd*Nws*Nwt,Nwd*Nws*Nwt);
for s=1:Nws
    sparsityAlongOneWindDir((s-1)*Nwt+1:s*Nwt, (s-1)*Nwt+1:s*Nwt) = tril(ones(Nwt,Nwt));
end
for w=1:Nwd
    sparsOnsetWindCnstrWrtU((w-1)*Nws*Nwt+1:w*Nws*Nwt, (w-1)*Nws*Nwt+1:w*Nws*Nwt) = sparsityAlongOneWindDir;
end

jacobianSparsityConstraints = [ zeros(1,2*Nwt), ones(1,Nwd*Nws*Nwt) ; ...
                    horzcat(relDjacobiSparsity,relDjacobiSparsity,zeros(numRelDcon,Nwd*Nws*Nwt));
                    horzcat(ones(Nwd*Nws*Nwt,2*Nwt),sparsOnsetWindCnstrWrtU) ];

end

%
% end findSparsityPattern.m
%