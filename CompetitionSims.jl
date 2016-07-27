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

    #Initialize parameters
    #k = 0.01*0.4*500.0*10^6; #Amazon estimated
    #k = 0.1*0.4*500.0*10^6; #Amazon klarge
    #k = 0.001*0.4*500.0*10^6; #Amazon ksmall
    k = 1.383*10^6; #estimated from 0.2 people /km^2

    a_hm = avec[i];
    c = cvec[j];

    a_mh = a_hm/c; #Effect of humans in Amazon

    T, den = ode45(frank_comp, initial_SA, [0., 10.0^7]);
    den = hcat(den...)';

    #Save trajectory
    T_am = T;
    den_am = den;

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

      #Re-initialize parameters
      a_mh = a_hm/c; #Humans do have an effect on monsters globally
      k = 10.0*10^9; #Global

      #Restart Global model with new variables
      initial = [1.01*10^9,Mout];
      T, den = ode45(frank_comp, initial, [Tout, 10.0^7]);
      den = hcat(den...)';

      T_glob = T;
      den_glob = den;

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
writedlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_text_lowdeath.csv",t_ext);

t_extam=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_textam.csv");
t_ext=readdlm("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/d_text.csv");

ext_plot=plot(
[layer(y=t_ext[:,j],x=avec, Geom.line, Theme(default_color=colorant"red")) for j in 1:lcvec]...,
[layer(y=t_extam[:,j],x=avec, Geom.line) for j in 1:lcvec]...,
Coord.Cartesian(xmin=1, xmax=10, ymin=1000,ymax=2*10^4),
Guide.xlabel("Monster aggressiveness"),Guide.ylabel("Time to extinction"))

draw(PDF("$(homedir())/Dropbox/PostDoc/2016_Frankenstein/fig_ext_AmCat.pdf", 8inch, 5inch), ext_plot)


#Combine T_am;;T_glob, den_am;;den_glob for a full trajectory
T_full = [T_am,T_glob];
den_full = [den_am,den_glob];
