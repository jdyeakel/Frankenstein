using ODE
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/src/frank_comp.jl")

a_hm = 5.5;
a_mh = a_hm/10;
k = 10.0*10^9;
r_h = 0.0067;
r_m = 0.0067;

initial = [1.09*10^9,2.0];

T, den = ode45(frank_comp, initial, [0., 10.0^4]);
den = hcat(den...)';

plot(
  layer(x=T,y=den[:,1],Geom.line),
  layer(x=T,y=den[:,2],Geom.line)
  )

########################

avec = collect(1.0:0.05:10.0);
lavec = length(avec);
cvec = collect(1.0:1.0:10.0);
lcvec = length(cvec);
t_ext = zeros(lavec,lcvec);

k = 10.0*10^9; #Global
#k = 500.0*10^6 #South America
r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
initial = [1.09*10^9,2.0]; #Global
#initial = [11.0*10^6,2.0]; #South America
a_hm=0.0;
a_mh=0.0;
c=0.0;
for j=1:lcvec
  for i=1:lavec
    #Initialize parameters
    a_hm = avec[i];
    c = cvec[j];

    a_mh = a_hm/c;

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
  print(j)
end





#
# ext_plot=plot([layer(y=t_ext[:,j],x=avec, Geom.line) for j in 1:lcvec]...,
# Coord.Cartesian(xmin=1, xmax=10, ymin=5000,ymax=1*10^5))
#
# draw(PDF("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/fig_ext_SA.pdf", 8inch, 5inch), ext_plot)


############## The European catalyst

avec = collect(1.0:0.05:10.0);
lavec = length(avec);
cvec = collect(1.0:1.0:10.0);
lcvec = length(cvec);
t_extEU = zeros(lavec,lcvec);


r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
#initial = [1.09*10^9,2.0]; #Global
initial_EU = [178.0*10^6,2.0]; #Europe
a_hm=0.0;
a_mh=0.0;
c=0.0;
k=0.0;
for j=1:lcvec
  for i=1:lavec
    a_hm = avec[i];
    c = cvec[j];

    a_mh = a_hm/c; #Effect of humans in Amazon


    #Initialize parameters
    k = 600*10^6; #from modern pop trajectory

    T, den = ode45(frank_comp, initial_EU, [0., 10.0^7]);
    den = hcat(den...)';


    hden=den[:,1];
    mden=den[:,2];
    #When do monsters outnumber humans? At this point, we assume they will expand...
    test_func = mden-hden;
    pos_outnum = find(x->x>0,test_func);
    if length(pos_outnum)==0
      t_extEU[i,j] = Inf;
    else
      Tout=T[pos_outnum[1]];
      Mout=mden[pos_outnum[1]];
      Hout=hden[pos_outnum[1]];
      #Hloss = initial_SA[1]-Hout;

      #Save trajectory
      T_EU = copy(T[1:pos_outnum[1]]);
      den_EU = copy( hcat(den[1:pos_outnum[1],1],den[1:pos_outnum[1],2]));

      #Re-initialize parameters
      a_mh = a_hm/c; #Humans do have an effect on monsters globally
      k = 10.0*10^9; #Global

      #Restart Global model with new variables
      #initial = [1.01*10^9,copy(Mout)];
      initial = [k,copy(Mout)]; #Start at Global carrying capacity

      T, den = ode45(frank_comp, initial, [Tout, 10.0^7]);
      den = hcat(den...)';

      T_glob = copy(T);
      den_glob = copy(den);

      #When does the human trajectory fall below threshold value?
      hden = den[:,1];
      pos = find(x->x<thresh,hden);
      if length(pos) > 0
        t_extEU[i,j] = T[pos[1]];
      else
        t_extEU[i,j] = Inf;
      end
    end
  end #end i
  print(j)
end #end j


writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_textEU_lowdeath.csv",t_extEU);




############## The south american catalyst
#SA is 6.9 million sq miles
#1816 density is 1.5 people per square mile
#2016 density is 5.8 people per square mile
#Amazon basin is 2.67 million sq mile (40% of total)
#say Amazon k = 0.4*500*10^6
#An evenly distributed populalation is 0.4*11*10^6...
#Let's say the population is only 5% of this 0.01*0.4*11*10^6
#1) explode in amazon basin
#2) world simulation
########################

avec = collect(1.0:0.05:10.0);
lavec = length(avec);
cvec = collect(1.0:1.0:10.0);
lcvec = length(cvec);
t_extam = zeros(lavec,lcvec);


r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
#initial = [1.09*10^9,2.0]; #Global
initial_SA = [(0.01*0.4*11.0*10^6),2.0]; #South America
a_hm=0.0;
a_mh=0.0;
c=0.0;
k=0.0;
for j=1:lcvec
  for i=1:lavec

    a_hm = avec[i];
    c = cvec[j];

    a_mh = a_hm/c; #Effect of humans in Amazon

    #Initialize parameters
    #k = 0.01*0.4*500.0*10^6; #Amazon estimated
    #k = 0.1*0.4*500.0*10^6; #Amazon klarge
    #k = 0.001*0.4*500.0*10^6; #Amazon ksmall
    k = 1.383*10^6; #estimated from 0.2 people /km^2

    T, den = ode45(frank_comp, initial_SA, [0., 10.0^7]);
    den = hcat(den...)';

    hden=den[:,1];
    mden=den[:,2];
    #When do monsters outnumber humans? At this point, we assume they will expand...
    test_func = mden-hden;
    pos_outnum = find(x->x>0,test_func);
    if length(pos_outnum)==0
      t_extam[i,j] = Inf;
    else
      Tout=T[pos_outnum[1]];
      Mout=mden[pos_outnum[1]];
      Hout=hden[pos_outnum[1]];
      #Hloss = initial_SA[1]-Hout;

      #Save trajectory
      T_am = copy(T[1:pos_outnum[1]]);
      den_am = copy( hcat(den[1:pos_outnum[1],1],den[1:pos_outnum[1],2]));

      #Re-initialize parameters
      a_mh = a_hm/c; #Humans do have an effect on monsters globally
      k = 10.0*10^9; #Global

      #Restart Global model with new variables
      #initial = [1.01*10^9,copy(Mout)];
      initial = [k,copy(Mout)]; #Start at Global carrying capacity

      T, den = ode45(frank_comp, initial, [Tout, 10.0^7]);
      den = hcat(den...)';

      T_glob = copy(T);
      den_glob = copy(den);

      #When does the human trajectory fall below threshold value?
      hden = den[:,1];
      pos = find(x->x<thresh,hden);
      if length(pos) > 0
        t_extam[i,j] = T[pos[1]];
      else
        t_extam[i,j] = Inf;
      end
    end
  end #end i
  print(j)
end #end j

writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_textam2_lowdeath.csv",t_extam);


t_extam=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_textam2_lowdeath.csv");
t_ext=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_textEU_lowdeath.csv");

ext_plot=plot(
[layer(y=t_ext[:,j],x=avec, Geom.line, Theme(default_color=colorant"red")) for j in 1:lcvec]...,
[layer(y=t_extam[:,j],x=avec, Geom.line, Theme(default_color=colorant"blue")) for j in 1:lcvec]...,
Coord.Cartesian(xmin=1, xmax=10, ymin=1000,ymax=1.5*10^4),
Guide.xlabel("Monster aggressiveness"),Guide.ylabel("Time to extinction"))

draw(PDF("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/fig_ext_AmCat.pdf", 5inch, 5inch), ext_plot)








#Full trajectory analysis
#For both European and Amazon catalyst scenarios

#European Catalyst
r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
#initial = [1.09*10^9,2.0]; #Global

initial_EU = [178.0*10^6,2.0]; #Europe

a_hm=8.0;
c = 10.0;

a_mh = a_hm/c; #Effect of humans in EU


#Initialize parameters
k = 600*10^6; #from modern pop trajectory

T, den = ode45(frank_comp, initial_EU, [0., 10.0^7]);
den = hcat(den...)';


hden=den[:,1];
mden=den[:,2];
#When do monsters outnumber humans? At this point, we assume they will expand...
test_func = mden-hden;
pos_outnum = find(x->x>0,test_func);
if length(pos_outnum)==0
  t_extEU[i,j] = Inf;
else
  Tout=T[pos_outnum[1]];
  Mout=mden[pos_outnum[1]];
  Hout=hden[pos_outnum[1]];
  #Hloss = initial_SA[1]-Hout;

  #Save trajectory
  T_EU = copy(T[1:pos_outnum[1]]);
  den_EU = copy( hcat(den[1:pos_outnum[1],1],den[1:pos_outnum[1],2]));

  #Re-initialize parameters
  a_mh = a_hm/c; #Humans do have an effect on monsters globally
  k = 10.0*10^9; #Global

  #Restart Global model with new variables
  #initial = [1.01*10^9,copy(Mout)]; #Starting at 1816 global densities
  initial = [k,copy(Mout)]; #Starting at Global carrying capacity

  T, den = ode45(frank_comp, initial, [Tout, 10.0^7]);
  den = hcat(den...)';

  T_glob = copy(T);
  den_glob = copy(den);

#Combine T_am;;T_glob, den_am;;den_glob for a full trajectory
T_full = [T_EU;T_glob];
den_full = [den_EU;den_glob];

TdenEU = hcat(T_full,den_full);
end
writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/traj_EUa2c4.csv",TdenEU);



#Amazon Catalyst
r_h = 0.0067;
r_m = 0.0067*1.5;
thresh = 1000.0;
#initial = [1.09*10^9,2.0]; #Global

initial_SA = [(0.01*0.4*11.0*10^6),2.0]; #South America


# a_hm=3.0;
# c = 4.0;
#
# a_mh = a_hm/c; #Effect of humans in Amazon

#Initialize parameters
#k = 0.01*0.4*500.0*10^6; #Amazon estimated
#k = 0.1*0.4*500.0*10^6; #Amazon klarge
#k = 0.001*0.4*500.0*10^6; #Amazon ksmall
k = 1.383*10^6; #estimated from 0.2 people /km^2

T, den = ode45(frank_comp, initial_SA, [0., 10.0^7]);
den = hcat(den...)';

hden=den[:,1];
mden=den[:,2];
#When do monsters outnumber humans? At this point, we assume they will expand...
test_func = mden-hden;
pos_outnum = find(x->x>0,test_func);
if length(pos_outnum)==0
  t_extam[i,j] = Inf;
else
  Tout=T[pos_outnum[1]];
  Mout=mden[pos_outnum[1]];
  Hout=hden[pos_outnum[1]];
  #Hloss = initial_SA[1]-Hout;

  #Save trajectory
  T_am = copy(T[1:pos_outnum[1]]);
  den_am = copy( hcat(den[1:pos_outnum[1],1],den[1:pos_outnum[1],2]));

  #Re-initialize parameters
  a_mh = a_hm/c; #Humans do have an effect on monsters globally
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

TdenAM = hcat(T_full,den_full);
end
writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/traj_ama2c4.csv",TdenAM);


TdenAM=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/traj_ama2c4.csv");
TdenEU=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/traj_EUa2c4.csv");



traj_plot = plot(
#Amazon Catalyst
layer(x=TdenAM[:,1],y=TdenAM[:,2],Geom.point,Theme(default_color=colorant"blue",default_point_size=1.4pt, highlight_width = 0pt)),
layer(x=TdenAM[:,1],y=TdenAM[:,2],Geom.line,Theme(default_color=colorant"blue")),

layer(x=TdenAM[:,1],y=TdenAM[:,3],Geom.point,Theme(default_color=colorant"green",default_point_size=1.4pt, highlight_width = 0pt)),
layer(x=TdenAM[:,1],y=TdenAM[:,3],Geom.line,Theme(default_color=colorant"green")),

#European Catalyst
layer(x=TdenEU[:,1],y=TdenEU[:,2],Geom.line,Theme(default_color=colorant"blue")),
layer(x=TdenEU[:,1],y=TdenEU[:,3],Geom.line,Theme(default_color=colorant"green")),
Scale.x_log10,Scale.y_log10,Coord.cartesian(ymin=0,ymax=12,xmin=2,xmax=4.0));
draw(PDF("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/fig_trajEUAM.pdf", 5inch, 5inch), traj_plot)
