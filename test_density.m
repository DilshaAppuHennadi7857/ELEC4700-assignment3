close all
clearvars

%%
%
% An example of finding the density of particles along the x dimension
%

% x defines the particles' positions along the x dimension
x = [80,10,60,70];

% nx defines the number of boxes the area is divided into (num bins if you
% will)
nx = 4;

% L is the actual length of the region
L = 200;

% If the length of the region is 200 and we divide it into 4 boxes, the bins
% will go as follows:
% x: 0   50  100 150 200
% ni:1   2   3   4 --> indexing of the bins

% we can determine the index that each particle falls into by dividing it's
% current position by the length of the area and multiplying by the number
% of bins
% NOTE: Be wary of the 0 case (it will yield a bin number of 0 which is
% invalid), we want to check for the 0 case and be sure those particles are
% placed in bin 1.
xi = ceil(x/L*nx)

% We will then go through the bins in the x direction and check to see
% which particles are in each bin. This is done through a logic array which
% gives us the sum of particles located in that bin.
% NOTE: This is only a 1D example, for the assignment, we are working in
% 2D. So we will have:
% a nested for loop for nx and ny
% matchx will use an AND operation to check for the x and y index
% ebox will be a matrix, and indexed using i,j
for i = 1:nx    
    matchx = xi==i;
    sum_e = sum(matchx);
    ebox(i) = sum_e;
end
ebox
A = [ 2 0 0 0; 0 1 0 0; 0 0 0 0]

% Once the matrix ebox is created, the density plot can just be plotted as
% a surface plot
surf(A)