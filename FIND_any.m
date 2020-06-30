function [R,C]=FIND_any(DM,sizex,sizey)

%Find a position where there is cell of any type on domain DM

R=randsample(sizex,1);
C=randsample(sizey,1);
while DM(R,C)==0
    R=randsample(sizex,1);
    C=randsample(sizey,1);
end