function [R,C]=FIND(x,DM,sizex,sizey)

%Find a position where there is cell of type x on domain DM

R=randsample(sizex,1);
C=randsample(sizey,1);
while DM(R,C)~=x
    R=randsample(sizex,1);
    C=randsample(sizey,1);
end