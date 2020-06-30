function [a_m, r_m, a_l, r_l, a_x, r_x, a_i, r_i, rand_move]= attraction_rates(mutant_type)

%This function stores the attraction rates for each of the different
%cell-cell interactions

%Outputs:
%a_m: vector of weight of attraction of m towards other m (a_mm), x (a_mx), id,(a_mi)
%r_m: vector of weight of repulsion of m towards other m (r_mm), x (r_mx), id,(r_mi)
%a_l: vector of weight of attraction of il towards m (a_lm), x (a_lx), id,(a_li)
%r_l: vector of weight of repulsion of il towards m (r_lm), x (r_lx), id,(r_li)
%a_x: vector of weight of attraction of x towards m (a_xm), x (a_xx), id,(a_xi)
%r_x: vector of weight of repulsion of x towards m (r_xm), x (r_xx), id,(r_xi)
%a_i: vector of weight of attraction of id towards m (a_im), x (a_ix), id,(a_ii)
%r_i: vector of weight of repulsion of id towards m (r_im), x (r_ix), id,(r_ii)
%rand_move : vector for a random move

%%RATES
a_mm=0;
a_mx=0;
a_mi=0;
r_mm=0;
r_mx=100;
r_mi=30;



a_xm=10; 
a_xx=0;
a_xi=1000;
r_xm=5; 
r_xx=0; 
r_xi=0;
a_im=10;
a_ix=300;
a_ii=0;
r_im=10;
r_ix=0;
r_ii=0;

a_lm=50;
a_lx=0;
a_li=50;
r_lm=50;
r_lx=50;
r_li=100;


if strcmp(mutant_type,'leo_nac')==1 || strcmp(mutant_type, 'leo')==1 ||  strcmp(mutant_type, 'leo_pfe')==1 || strcmp(mutant_type, 'leo_shd')...
        || strcmp(mutant_type, 'leo_sbr')==1
    r_mx=0;
end


a_m=[a_mm;a_mx;a_mi];
r_m=[r_mm;r_mx;r_mi];


r_l=[r_lm;r_lx;r_li];
a_l=[a_lm;a_lx;a_li];

a_x=[a_xm;a_xx;a_xi];
r_x=[r_xm;r_xx;r_xi];

a_i=[a_im;a_ix;a_ii];
r_i=[r_im;r_ix;r_ii];

Pb=1;
rand_move=[1;1;1;1;1;1;1;1];
end






