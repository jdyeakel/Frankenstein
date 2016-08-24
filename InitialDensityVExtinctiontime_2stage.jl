using ODE
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")

#Simulate extinction time across different initial human population sizes

den_power = collect(0.5:0.01:10);
l_den = length(den_power);
t_ext = zeros(l_den);

r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
a_hm=2.0;
c = 3.0;
a_mh = a_hm/c; #Effect of humans
#Initialize carrying capacity carrying capacity
k = 0.0
for i in 1:l_den
  power = den_power[i];
  initial_density = [10^power,2.0]; #Europe
  # Here, the last iteration equals k

  k = 1.383*10^6; #estimated from 0.2 people /km^2

  T, den = ode45(frank_comp, initial_density, [0., 10.0^7]);
  den = hcat(den...)';

  hden=den[:,1];
  mden=den[:,2];
  #When do monsters outnumber humans? At this point, we assume they will expand...
  test_func = mden-(0.9*k);
  pos_outnum = find(x->x>=0,test_func);
  if length(pos_outnum)==0
    t_ext[i] = Inf;
  else
    Tout=copy(T[pos_outnum[1]]);
    Mout=copy(mden[pos_outnum[1]]);
    Hout=copy(hden[pos_outnum[1]]);
    #Hloss = initial_SA[1]-Hout;

    #Save trajectory
    T_am = copy(T[1:pos_outnum[1]]);
    den_am = copy( hcat(den[1:pos_outnum[1],1],den[1:pos_outnum[1],2]));

    #The new global carrying capacity
    k = 10.0*10^9; #Global

    #Restart Global model with new variables
    #initial = [1.01*10^9,copy(Mout)]; #start at Global 1816 density
    initial = [k,copy(Mout)]; #Start at Global carrying capacity

    T, den = ode45(frank_comp, initial, [Tout, 10.0^7]);
    den = hcat(den...)';

    T_glob = copy(T);
    den_glob = copy(den);


    #Combine T_am;;T_glob, den_am;;den_glob for a full trajectory
    T_full = [T_am;T_glob];
    den_full = [den_am;den_glob];


    #When does the human trajectory fall below threshold value?
    hden = den_glob[:,1];
    pos = find(x->x<thresh,hden);
    if length(pos) > 0
      t_ext[i] = T_glob[pos[1]];
    else
      t_ext[i] = Inf;
    end
  end
  print(i);
end

writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/initalvext2stage.csv",t_ext);

den_actual = [10^den_power[i] for i=1:l_den];
initial_ext_plot=plot(x=den_actual,y=t_ext,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt),Scale.x_log10,Scale.y_log10,
Coord.cartesian(xmin=2,xmax=6,ymin=5.8,ymax=7),Guide.xlabel("Initial human density"),Guide.ylabel("Time to extinction"));

initial_ext_plot=plot(x=den_actual,y=t_ext,Scale.x_log10,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt),
Coord.cartesian(xmin=0,xmax=10,ymin=6*10^3,ymax=1.2*10^4),Guide.xlabel("Initial human density"),Guide.ylabel("Time to extinction"));
#

draw(PDF("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/fig_InitialvExtinct2stage.pdf", 3inch, 3inch), initial_ext_plot)

#ymin=0,ymax=6,



plot(
layer(x=T,y=den[:,1],Geom.line),
layer(x=T,y=den[:,2],Geom.line),
Coord.cartesian(ymin=10^1.0,ymax=10^10.0),Guide.xlabel("10^i"),Guide.ylabel("Time to extinction"))
