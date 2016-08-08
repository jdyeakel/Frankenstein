using ODE
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")

#Simulate extinction time across different initial human population sizes

den_power = collect([0:0.01:6]);
l_den = length(den_power);
t_ext = zeros(l_den);
for i in 1:length(den_power)
  power = den_power[i];
  initial_density = [10^power,2.0]; #Europe

  r_h = 0.0067;
  r_m = 0.0067*1.5;
  thresh = 1000.0;
  a_hm=8.0;
  c = 10.0;
  a_mh = a_hm/c; #Effect of humans in EU

  #Global carrying capacity
  k = 10.0*10^9; #Global

  T, den = ode45(frank_comp, initial_density, [Tout, 10.0^7]);
  den = hcat(den...)';

  #When does the human trajectory fall below threshold value?
  hden = den[:,1];
  pos = find(x->x<thresh,hden);
  if length(pos) > 0
    t_ext[i] = T[pos[1]];
  else
    t_ext[i] = Inf;
  end

end
