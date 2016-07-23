using ODE
using Gadfly
include("/home/jdyeakel/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")

a_hm = 1.5;
a_mh = a_hm/2;
k = 10*10^9;
r_h = 0.0067;
r_m = 0.0067;

initial = [1.09*10^9,1.0];

T, den = ode45(frank_comp, initial, [0., 10.0^7]);
den = hcat(den...)';

plot(
  layer(x=T,y=den[:,1],Geom.line),
  layer(x=T,y=den[:,2],Geom.line),
  Scale.y_log10
  )

########################

avec = collect(1.0:0.1:10.0);
lavec = length(avec);
cvec = collect(1.0:1:10.0)
lcvec = length(cvec);
t_ext = zeros(lavec,lcvec);
for j=1:lcvec
  for i=1:lavec
    #Initialize parameters
    a_hm = avec[i];
    c = cvec[j];
    a_mh = a_hm/c;
    k = 10*10^9;
    r_h = 0.0067;
    r_m = 0.0067;
    thresh = 1000.0;
    initial = [1.09*10^9,1.0];

    T, den = ode45(frank_comp, initial, [0., 10.0^7]);
    den = hcat(den...)';

    #When does the human trajectory fall below threshold value?
    hden = den[:,1];
    pos = find(x->x<thresh,hden);
    if length(pos) > 0
      t_ext[i,j] = T[pos[1]];
    else
      t_ext[i,j] = Inf;
    end
  end
end


plot([layer(y=t_ext[:,j],x=avec, Geom.line) for j in 1:lcvec]...,Coord.Cartesian(xmin=1, xmax=10))
