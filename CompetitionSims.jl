using ODE
using Gadfly

a_hm = 1.5;
a_mh = a_hm/2;
k = 10*10^9;
r_h = 0.0067;
r_m = 0.0067;

initial = [1.09*10^9,1.0];

include("/home/jdyeakel/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")
T, den = ode45(frank_comp, initial, [0., 10.0^7]);
den = hcat(den...)';

plot(
  layer(x=T,y=den[:,1],Geom.line),
  layer(x=T,y=den[:,2],Geom.line),
  Scale.y_log10
  )

for i=1:
