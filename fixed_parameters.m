function [short_distance_m, short_distance_x, max_distance, long_x, growth_rateX,...
    growth_rateY, a_m, r_m, a_l, r_l, a_x, r_x, a_i, r_i, rand_move]=...
    fixed_parameters(mutant_type)

%THESE PARAMETERS DO NOT CHANGE WITH TIME
site_size_m=40;

%Distances used in the model
short_distance_m=1;
short_distance_x=2;
max_distance=2;
long_x=12;

%Rate of growth in x and y axis
if strcmp(mutant_type,'ablation')==0 && strcmp(mutant_type,'kondo_exp')==0
    growth_rateX=33/(site_size_m*60*24);  %30um per day = 30/sitesize sites per day
    growth_rateY=130/(site_size_m*60*24);  %130um per day = 130/sitesize sites per day
else
    growth_rateX=0;
    growth_rateY=0;
end

%Attraction rates
[a_m, r_m, a_l, r_l, a_x, r_x, a_i, r_i, rand_move]= attraction_rates(mutant_type);

    



