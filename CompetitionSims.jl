using ODE
using Gadfly

a_mh = 0.5;
a_hm = 0.5;
k_m = 100;
k_h = 100;
r_h = 0.5;
r_m = 0.5;

initial = [20.0,10.0];

include("/home/jdyeakel/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")
T, den = ode23(frank_comp, initial, [0., 40]);
den = hcat(den...)';

plot(x=T,y=den[:,1],Geom.line)
