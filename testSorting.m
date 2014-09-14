clear all;close all;
Nwt=5;
xPos = [100;80;0;1800;0];
yPos = [100;300;90;1800;1800];
turbineSites = [xPos, yPos];

[indexedTurbines,I] = sortrows(turbineSites,[-2,1]);
[dumDum,rI]=sortrows(I,1);

%{
By the indexing I and reverse indexing rI we have that:
  turbineSites(I,:) == indexedTurbines
  indexedTurbines(rI,:) == turbineSites
%}