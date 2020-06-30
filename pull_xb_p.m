function [domain_matrix]=pull_xb_p(domain_matrix,R,C,sizex,sizey)

done=0;
P_x=[R+5,R-5,R,R,...
    R+1,R-1,R+1, R-1, ...
    R+2, R-2, R+2, R-2, ...
    R-3, R+3, R-3, R+3, ...
    R-4, R+4, R-4, R+4, ...
    R-5, R+5, R-5, R+5, ...
    R+5, R-5, R+5, R-5];

P_y=[C,C,C+5,C-5,...
    C+5, C+5, C-5, C-5, ...
    C+5, C+5, C-5, C-5, ...
    C-4, C+4, C+4, C-4, ...
    C-3, C+3, C-3, C+3, ...
    C+1, C-1, C+1, C-1, ...
    C+2, C+2, C-2, C-2];

r=randsample(length(P_x),length(P_x)); %Shuffle the options
k=1;

while k<=length(P_x) && done==0 %Go through each option sequentially
    M=P_x(r(k));
    N=P_y(r(k));
    if N>sizey || N<0
        N=mod(N,sizey); %periodic BCs
    elseif N==0
        N=sizey; %periodic BCs
    end
    k=k+1;
    if  M>0 && M<= sizex && domain_matrix(ceil(M/2),ceil(N/2))~=0
        done=1; %meets criteria! it's on the board and theres a melanophore nearby.
    end
end %choose a nearby melanophore

if done==1 && domain_matrix(ceil(R/2),ceil(C/2))==0
    %move the melanophore
    domain_matrix(ceil(M/2),ceil(N/2))=0;
    domain_matrix(ceil(R/2),ceil(C/2))=1;
end