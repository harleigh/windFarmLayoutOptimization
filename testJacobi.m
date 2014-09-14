%Test the jacobian wrt relative distance constraints
clear all; close all;

Nwt=5;
farmWidth=1800;
showFigure=0;
layoutScheme = empiricalScheme(farmWidth,Nwt,showFigure);
sparsity = find(triu(ones(Nwt,Nwt),1));

turbineSites = [layoutScheme(1:Nwt) layoutScheme(Nwt+1:2*Nwt)];
XLOC = 1;  %enumeration: location of xPositions in turbineSites is col 1
YLOC = 2;  %enumeration: location of yPositions in turbineSites is col 2
indexedTurbines = sortrows(turbineSites,-YLOC);


dwD = calcDownWindDist(indexedTurbines(:,YLOC));
cwD = calcCrossWindDist(indexedTurbines(:,XLOC));
relDsquared = calcSquaredRelativeDist(cwD, dwD);

relDsquared(sparsity);

[sparsityRow, sparsityCol] = ind2sub([Nwt,Nwt],sparsity);
totRelSquaredDistConstraints = numel(sparsity);
dumDum = zeros(totRelSquaredDistConstraints,Nwt);

for i=1:totRelSquaredDistConstraints
    dumDum(i,sparsityRow(i))=1;
    dumDum(i,sparsityCol(i))=1;
end

jacobianSparsity = horzcat(dumDum,dumDum);



