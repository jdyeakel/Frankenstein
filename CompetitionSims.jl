Using ODE

a_mh = 0.5;
a_hm = 0.5;
k_m = 100;
k_h = 100;
r_h = 0.5;
r_m = 0.5;



T, den = ode23(frank_comp, initial, [0., 40]);
den = hcat(den...)';

