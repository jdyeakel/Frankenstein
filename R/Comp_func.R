Comp_func <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dH <- Rh*H*(1 - ((H + a_hm*M)/Kh))
    dM <- Rm*M*(1 - ((M + a_mh*H)/Km))
    
    #Return rate of change
    list(c(dH,dM))
  }) #end with statement
  
}