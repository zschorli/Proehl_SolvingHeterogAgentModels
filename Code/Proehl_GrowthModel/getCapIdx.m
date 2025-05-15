% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function determines the grid index of the given grid which is 
% closest to the input cap.
%__________________________________________________________________________

function [capIdx,capIdx1,capIdx2] = getCapIdx(cap,grid)

if isrow(cap)
   cap=cap';
end

capIdx1 = sum(repmat(grid,[length(cap) 1])-...
          repmat(cap,[1 length(grid)])<=0,2);
capIdx2 = capIdx1+1;
capIdx1 = (capIdx1>=length(grid))*(length(grid)-1)+...
          (capIdx1<length(grid)).*capIdx1;
capIdx1 = (capIdx1<1)*1+(capIdx1>=1).*capIdx1;      
capIdx2 = (capIdx2>length(grid))*length(grid)+...
          (capIdx2<=length(grid)).*capIdx2;

capIdx = (abs(grid(capIdx1)'-cap)<=...
         abs(grid(capIdx2)'-cap)).*capIdx1+...
         (abs(grid(capIdx2)'-cap)<...
         abs(grid(capIdx1)'-cap)).*capIdx2;

end