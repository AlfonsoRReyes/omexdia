## =============================================================================
## R-code to run OMEXDIA from the DLL
## Implementation by karline soetaert (karline.soetaert@nioz.nl)
## =============================================================================

  N     <- 100
  Grid  <- setup.grid.1D(x.up=0, dx.1 = 0.01, N = N, L = 50)
  Depth <- Grid$x.mid                  # depth of each box

# porosity gradients
  exp.profile <- function(x, y.0, y.inf, x.att = 1, L = 0)
           return(y.inf + (y.0 - y.inf) * exp(-pmax(0, x-L)/x.att))

  porGrid <- setup.prop.1D(func = exp.profile,
                         grid = Grid,
                         y.0 = 0.9, y.inf = 0.5,
                         x.att = 3)

## The parameters
## the default parameter values

  Parms <- c(
  ## organic matter dynamics  #
    MeanFlux = 20000/12*100/365,  # nmol/cm2/d - Carbon deposition: 20gC/m2/yr
    rFast    = 0.07            ,  #/day        - decay rate fast decay detritus
    rSlow    = 0.00001         ,  #/day        - decay rate slow decay detritus
    pFast    = 0.9             ,  #-           - fraction fast detritus in flux
    w        = 0.1/1000/365    ,  # cm/d       - advection rate
    NCrFdet  = 0.16            ,  # molN/molC  - NC ratio fast decay detritus
    NCrSdet  = 0.13            ,  # molN/molC  - NC ratio slow decay detritus

  ## oxygen and DIN dynamics  #

  ## Nutrient bottom water conditions
    bwO2            = 300      ,    #mmol/m3     Oxygen conc in bottom water
    bwNO3           = 10       ,    #mmol/m3
    bwNH3           = 1        ,    #mmol/m3
    bwODU           = 0        ,    #mmol/m3

  ## Bioturbation
    biot            = 1/365    ,    # cm2/d      - bioturbation coefficient
    mixL            = 5        ,    # cm         - depth of mixed layer

  ## Nutrient parameters
    NH3Ads          = 1.3      ,    #-           Adsorption coeff ammonium
    rnit            = 20.      ,    #/d          Max nitrification rate
    ksO2nitri       = 1.       ,    #umolO2/m3   half-sat O2 in nitrification
    rODUox          = 20.      ,    #/d          Max rate oxidation of ODU
    ksO2oduox       = 1.       ,    #mmolO2/m3   half-sat O2 in oxidation of ODU
    ksO2oxic        = 3.       ,    #mmolO2/m3   half-sat O2 in oxic mineralisation
    ksNO3denit      = 30.      ,    #mmolNO3/m3  half-sat NO3 in denitrification
    kinO2denit      = 1.       ,    #mmolO2/m3   half-sat O2 inhib denitrification
    kinNO3anox      = 1.       ,    #mmolNO3/m3  half-sat NO3 inhib anoxic degr
    kinO2anox       = 1.       ,    #mmolO2/m3   half-sat O2 inhib anoxic min

  ## Diffusion coefficients, temp = 10dgC
    #Temp            = 10                     ,   # temperature
    DispO2          = 0.955    +10*0.0386    ,  #cm2/d
    DispNO3         = 0.844992 +10*0.0336    ,
    DispNH3         = 0.84672  +10*0.0336    ,
    DispODU         = 0.8424   +10*0.0242    ,
    TOC0            = 0.5)

## SHINY application

OMEXDIAsteady <- function (pars= list(), D = 60) {
  DIA <-  OMEXDIAsolve (pars, D)
  plotDIA (DIA)
}
  
OMEXDIAsolve <- function (pars = list(), D = 60)  {

## check parameter inputs
  nms <- names(Parms)
  Parms[(namc <- names(pars))] <- pars
  if (length(noNms <- namc[!namc %in% nms]) > 0)
    warning("unknown names in pars: ", paste(noNms, collapse = ", "))
  
# Bioturbation profile
  Db <- setup.prop.1D(func = exp.profile,
                         grid = Grid,
                         y.0 = Parms[["biot"]], y.inf = 0.,
                         L = Parms[["mixL"]])$int
  parms <- Parms[-which(nms %in%c("biot", "mixL"))]

  Flux <- as.double(Parms["MeanFlux"])

  initpar <- c(parms, Grid$dx, Grid$dx.aux, porGrid$mid, porGrid$int,Db)

## solve the steady-state condition
  OC   <- rep(10,6*N)
  nout <- 1002
  outnames <- c("O2flux", "O2deepflux", "NO3flux", "NO3deepflux", "NH3flux",
      "NH3deepflux", "ODUflux", "ODUdeepflux", "Cflux", "partDenit", 
      "partAnoxic", "partOxic", 
      rep("Cprod",N), rep("Nprod",N), rep("TOC",N), 
      rep("Oxicmin",N), rep("Denitrific",N), rep("anoxicmin",N),
      rep("nitri",N), rep("oduox",N), rep("odudepo",N))
  
  ynames <- c("FDET","SDET","O2","NO3","NH3","ODU")

  DIA  <- steady.1D(y=as.double(OC),func="omexdiamod",initfunc="initomexdia",
                     names = ynames, initforc = "initforc", forcings = mean(Flux),  
                     initpar=initpar,nspec=6,
                     dllname="omexdia", 
                     nout=nout,outnames = outnames, positive=TRUE)

  DIA$Parms <- Parms
  return(DIA)
   