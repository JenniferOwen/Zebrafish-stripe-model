function [total]=calc_neighbours(R,C,x,max_distance,sizex,sizey,DM,add,inc)

%Uses the uniform norm to calculate number of neighbours of type x on the
%domain matrix at distance max_distance from space R,C.
%Periodic BCs in the y direction
%non-periodic BCs in the x direction

%Input
%R: location of cell in x direction
%C: location of cell in y direction
%x: indicator of cell type
%max_distance: largest distance tested
%sizex: total size of domain in number of lattice spaces in x direction
%sizey: total size of domain in number of lattice spaces in y direction
%DM: domain matrix
%add: weight added per occupied site
%inc: determines if the area for neighbours is just at boundary (inc = 0)
%      or whether it includes the whole area (inc = 1)

%Output
%total: number of neighbours of type x in boundary max_distance away from
%centre cell R,C on domain matrix DM.

total=0;

for i=-max_distance:max_distance
    for j=-max_distance:max_distance
        if inc==1 || (inc==0 && (abs(i)==max_distance) || (abs(j)==max_distance))
                R_check=R+i; %choose correct R and C to check
                C_check=C+j;
            if (R_check>sizex || R_check<1 ) || (R_check==R && C_check==C)
                continue %periodic BCs and do not include the site that the cell is already occupying.
            end
            if C_check>sizey || C_check<0
                C_check=mod(C_check,sizey); %periodic BCs
            elseif C_check==0
                C_check=sizey; %periodic BCs
            end
            if x==1
                if DM(R_check,C_check)~=0
                    total=total+add;
                end
            else
                if DM(R_check,C_check)==x
                    total=total+add;
                end
            end
        end
    end
    %     end
end



