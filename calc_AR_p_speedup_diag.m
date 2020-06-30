function [Attract,Repulse]=calc_AR_p_speedup_diag(R,C,x_num,m_num,irr_num,max_distance,sizex,sizey, domain_matrix, domain_matrix_X, domain_matrix_irr)

%Function checks the surroundings of R,C and outputs matrices of the
%numbers of cells in each position - to be used with the weighting vector
%(based on attraction/repulsion to cells to determine probability of  
%moving in each given direction

%Inputs:
%[R,C]: position of interest (POI)
%x_num: cell number of interest of the x domain (either 2 or 4 for xb or x)
%m_num: cell number of interest on the m domain (always 1)
%irr_num: cell number of interest on the i domain (either 5 or 6 for id or il)
%max_distance: distance of interest in number of lattice sites on domain of POI
%[sizex,sizey]: size of the domain where POI is located
%domain_matrix: M domain
%domain_matrix_X: X domain
%domain_matrix_irr: I domain 

%Outputs (to be used with the weighting matrix):
%Attract: number of cells of M, X, I chosen of interest around position (R,C)
%Repulse: number of cells of M, X, I chosen of interest in the opposite
%direction

%non-periodic BCs in the y direction
%periodic BCs in the x direction

%Counts the numbers of melanocytes in each direction
Nm=0;
Sm=0;
Em=0;
Wm=0;
SWm=0;
NWm=0;
SEm=0;
NEm=0;

%Counts the numbers of xanthophores in each direction
Nx=0;
Sx=0;
Ex=0;
Wx=0;
SWx=0;
NWx=0;
SEx=0;
NEx=0;

%Counts the numbers of dense iridophores in each direction
Ni=0;
Si=0;
Ei=0;
Wi=0;
SWi=0;
NWi=0;
SEi=0;
NEi=0;

for i=-max_distance:1:max_distance
    for j=-max_distance:1:max_distance
        R_check=R+i; %choose correct R and C to check
        C_check=C+j;
        if R_check>sizex || R_check<1 
            continue %non-periodic BCs
        end
        if C_check>sizey || C_check<0
            C_check=mod(C_check,sizey); %periodic BCs
        elseif C_check==0
            C_check=sizey; %periodic BCs
        end
        
        store_x= domain_matrix_X(R_check,C_check)==x_num;
        store_i= domain_matrix_irr(R_check,C_check)==irr_num;
        store_m= domain_matrix(ceil(R_check/2), ceil(C_check/2))~=0;
        if i<=0
            if i<0
                Nm=Nm+store_m;
                Nx=Nx+store_x;
                Ni=Ni+store_i;
            end
            
            if j>=0 && i+j~=0
                NEm=NEm+store_m;
                NEx=NEx+store_x;
                NEi=NEi+store_i;
            end
            
            if j<=0 && i+j~=0
                NWm=NWm+4;
                NWx=NWx+store_x;
                NWi=NWi+store_i;
            end
        end
        if i>=0
            if i>0
                Sm=Sm+store_m;
                Sx=Sx+store_x;
                Si=Si+store_i;
            end
            
            if j>=0 && i+j~=0
                SEm=SEm+store_m;
                SEx=SEx+store_x;
                SEi=SEi+store_i;
            end
            
            if j<=0 && i+j~=0
                SWm=SWm+store_m;
                SWx=SWx+store_x;
                SWi=SWi+store_i;
            end
        end
        if j<0
            Wm=Wm+store_m;
            Wx=Wx+store_x;
            Wi=Wi+store_i;
        end
        if j>0
            Em=Em+store_m;
            Ex=Ex+store_x;
            Ei=Ei+store_i;
        end
    end
    
end
Attract=[Em,Ex,Ei; Wm,Wx,Wi; Nm,Nx,Ni; Sm,Sx,Si; NEm,NEx,NEi; NWm,NWx,NWi; SEm,SEx,SEi; SWm,SWx,SWi];

Repulse=[Wm,Wx,Wi; Em,Ex,Ei; Sm,Sx,Si; Nm,Nx,Ni; SWm,SWx,SWi; SEm,SEx,SEi; NWm,NWx,NWi; NEm,NEx,NEi;];
end




