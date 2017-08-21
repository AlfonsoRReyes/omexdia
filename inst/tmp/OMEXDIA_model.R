## =============================================================================
##
## OMEXDIA_C: C, N, and O2 diagenesis
##
## The model function, parameters, grids, budget functions 
##
## Implementation by karline soetaert (karline.soetaert@nioz.nl)
## =============================================================================

#====================#
# Model equations    #
#====================#

OMEXDIA_C <-function(time = 0,     # time 
                     state,        # FDET, SDET, O2, NO3, NH3, ODU
                     parms,        # parameter values
                     ...)  {       

with (as.list(parms), {

  FDET  <- state[1:N]
  SDET  <- state[(N+1)  :(2*N)]
  O2    <- state[(2*N+1):(3*N)]
  NO3   <- state[(3*N+1):(4*N)]
  NH3   <- state[(4*N+1):(5*N)]
  ODU   <- state[(5*N+1):(6*N)]

# ----------------------------------------
# Transport
# ----------------------------------------   
# C deposition
# For steady-state, this will be = MeanFlux, for transient simulation, a sine
  Flux  <- MeanFlux * (1+sin(2*pi*time/365))     

# Solid substances
  FDETtran <- tran.1D(C = FDET, flux.up = Flux*pFast    , D = DbGrid, v = w,
                     AFDW = 0.5, VF = porGridSolid, dx = Grid)

  SDETtran <- tran.1D(C = SDET, flux.up = Flux*(1-pFast), D = DbGrid, v = w,
                     AFDW = 0.5, VF = porGridSolid, dx = Grid)

# Dissolved substances
  O2tran  <- tran.1D (C = O2, C.up = bwO2, D = DispO2, v = w,
                     VF = porGrid, dx = Grid)

  NO3tran <- tran.1D (C = NO3, C.up = bwNO3, D = DispNO3, v = w,
                     VF = porGrid, dx = Grid)

  NH3tran <- tran.1D (C = NH3, C.up = bwNH3, D = DispNH3/(1+NH3Ads),
                     v = w, VF = porGrid, dx = Grid)

  ODUtran <- tran.1D (C = ODU, C.up = bwODU, D = DispODU, v = w,
                     VF = porGrid, dx = Grid)

# ----------------------------------------
# Biogeochemical transformations
# ----------------------------------------

# porosity factors
  p2liquid    <- (1.-porGrid$mid)/porGrid$mid        # from solid -> liquid
  p2solid     <- porGrid$mid/(1.-porGrid$mid)        # from liquid -> solid  

# production of DIC, DIN, expressed per cm3 LIQUID/day
  DICprod_Min <- (rFast*FDET         + rSlow*SDET       ) * p2liquid
  DINprod_Min <- (rFast*FDET*NCrFdet + rSlow*SDET*NCrSdet)* p2liquid

# oxic mineralisation, denitrification, anoxic mineralisation
  Oxicminlim <- O2/(O2 + ksO2oxic)                   # limitation terms
  Denitrilim <- (1 - O2/(O2+kinO2denit))* NO3/(NO3+ksNO3denit)
  Anoxiclim  <- (1 - O2/(O2+kinO2anox ))*(1 - NO3/(NO3+kinNO3anox))
  Rescale    <- 1/(Oxicminlim+Denitrilim+Anoxiclim)

  OxicMin    <- DICprod_Min * Oxicminlim * Rescale   # oxic mineralisation
  Denitrific <- DICprod_Min * Denitrilim * Rescale   # Denitrification
  AnoxicMin  <- DICprod_Min * Anoxiclim  * Rescale   # an0xic mineralisation

# reoxidation and ODU deposition
  Nitri      <- rnit   * NH3 * O2/(O2+ksO2nitri)
  OduOx      <- rODUox * ODU * O2/(O2+ksO2oduox)

  pDepo      <- min(1,0.233*(w*365)^0.336 )
  OduDepo    <- AnoxicMin*pDepo

# ----------------------------------------
# Update the rate of change 
# ----------------------------------------

  dFDET <- FDETtran$dC - rFast*FDET  
  dSDET <- SDETtran$dC - rSlow*SDET

  dO2   <- O2tran$dC  -  OxicMin      -2* Nitri -      OduOx 

  dNH3  <- NH3tran$dC + (DINprod_Min  - Nitri ) / (1.+NH3Ads) 

  dNO3  <- NO3tran$dC - 0.8*Denitrific + Nitri

  dODU  <- ODUtran$dC + AnoxicMin  - OduOx - OduDepo

# ----------------------------------------
# output variables
# ----------------------------------------

  Norgflux      <- Flux*pFast*NCrFdet + Flux*(1-pFast)*NCrSdet
  Norgdeepflux  <- FDETtran$flux.down*NCrFdet + SDETtran$flux.down*NCrSdet
  NH3adsorption <- (DINprod_Min  - Nitri) * (1-1/(1+NH3Ads))
  TOC           <- (FDET+SDET)*1200/10^9/2.5  # % organic carbon (excess)

 # Model output
  return(list(
    c(dFDET, dSDET, dO2, dNO3, dNH3, dODU), # derivatives
    TOC  = TOC,                       # % organic carbon (excess)
 
 # exchange fluxes
    Norgflux = Norgflux,
    Norgdeepflux = Norgdeepflux,
    O2flux = O2tran$flux.up,          # O2 sediment-water flux
    O2deepflux = O2tran$flux.down,    # O2 deep(burial) flux
    NO3flux = NO3tran$flux.up,        # NO3 sediment-water flux
    NO3deepflux = NO3tran$flux.down,  # NO3 deep(burial) flux
    NH3flux = NH3tran$flux.up*(1+NH3Ads),        # NH3 sediment-water flux 
    NH3deepflux = NH3tran$flux.down*(1+NH3Ads),  # NH3 deep (burial) flux
    ODUflux = ODUtran$flux.up,        # ODU sediment-water flux
    ODUdeepflux = ODUtran$flux.down,  # ODU deep(burial) flux

 # rate profiles
    NH3adsorption = NH3adsorption,
    OxicMin = OxicMin,                # oxic mineralisation
    Denitrific = Denitrific,          # denitrification rates
    AnoxicMin = AnoxicMin,            # anoxic mineralisation
    Nitri = Nitri,                    # nitrification rates
    OduOx = OduOx))                   # ODU oxidation rates

})

} # end of OMEXDIA_C

require(ReacTran)   # The package with the solution methods
require(marelac)    # Tool for aquatic sciences

# reference abiotic conditions
Temp <- 9.5
Sal  <- 29
por  <- 0.4

#====================#
# Model applications #
#====================#

# Grid: 100 layers; total length=50 cm, first box=0.01 cm
Grid  <- setup.grid.1D(N = 100, dx.1 = 0.01, L = 50)
Depth <- Grid$x.mid
N     <- Grid$N

# porosity: fixed and constant value
porGrid <- setup.prop.1D(value = por, grid = Grid)

# 1 - porosity
porGridSolid <- setup.prop.1D(value = 1 - por, grid = Grid )

# Diffusion coefficients, uses diffcoeffs from package marelac
DiffCoeffs <- diffcoeff(S = Sal, t = Temp)*3600*24*1e4 # from m2/s -> cm2/d
DispO2     <- as.numeric(DiffCoeffs["O2"]  ) * por
DispNO3    <- as.numeric(DiffCoeffs["NO3"] ) * por
DispNH3    <- as.numeric(DiffCoeffs["NH3"] ) * por
DispODU    <- as.numeric(DiffCoeffs["HSO4"]) * por

# exponential function
exp.profile <- function(x, x.0, y.0, y.inf, x.att)
     return(y.inf + (y.0-y.inf)*exp(-pmax(0.,(x-x.0))/x.att))

# Bioturbation
biot <- 1/365           # cm2/d      - bioturbation coefficient
mixL <- 10              # cm         - depth of mixed layer

DbGrid <- setup.prop.1D(func = exp.profile, x.0 = mixL,
                      y.0 = biot, y.inf = 0, x.att = 1, 
                      grid = Grid)

# organic matter dynamics  #
MeanFlux <- 1000/12*100/365  # nmol/cm2/d - C deposition: 5gC/m2/yr
rFast    <- 0.1              #/day        - decay rate fast decay det.
rSlow    <- 0.01             #/day        - decay rate slow decay det.
pFast    <- 0.5              #-           - fraction fast det. in flux
w        <- 0.5/365          # cm/d       - advection rate
NCrFdet  <- 0.156            # molN/molC  - NC ratio fast decay det.
NCrSdet  <- 0.156            # molN/molC  - NC ratio slow decay det.

# Nutrient bottom water conditions
bwO2       <- 300       #mmol/m3     Oxygen conc in bottom water
bwNO3      <- 16        #mmol/m3
bwNH3      <- 7.5       #mmol/m3
bwODU      <- 0         #mmol/m3   

# Nutrient parameters
NH3Ads     <- 1.3    #-           Adsorption coeff ammonium
rnit       <- 2.     #/d          Max nitrification rate
ksO2nitri  <- 1.     #umolO2/m3   half-sat O2 in nitrification
rODUox     <- 2.     #/d          Max rate oxidation of ODU
ksO2oduox  <- 1.     #mmolO2/m3   half-sat O2 in oxidation of ODU
ksO2oxic   <- 3.     #mmolO2/m3   half-sat O2 in oxic mineralis
ksNO3denit <- 30.    #mmolNO3/m3  half-sat NO3 in denitrif
kinO2denit <- 1.     #mmolO2/m3   half-sat O2 inhib denitrif
kinNO3anox <- 1.     #mmolNO3/m3  half-sat NO3 inhib anoxic min
kinO2anox  <- 1.     #mmolO2/m3   half-sat O2 inhib anoxic min

pars <- c(
  MeanFlux   =  MeanFlux   ,rFast     =  rFast      ,
  rSlow      =  rSlow      ,pFast     =  pFast      ,
  w          =  w          ,NCrFdet   =  NCrFdet    ,
  NCrSdet    =  NCrSdet    ,bwO2      =  bwO2       ,
  bwNO3      =  bwNO3      ,bwNH3     =  bwNH3      ,
  bwODU      =  bwODU      ,NH3Ads    =  NH3Ads      ,
  rnit       =  rnit       ,ksO2nitri =  ksO2nitri  ,
  rODUox     =  rODUox     ,ksO2oduox =  ksO2oduox  ,
  ksO2oxic   =  ksO2oxic   ,ksNO3denit=  ksNO3denit ,
  kinO2denit =  kinO2denit ,kinNO3anox=  kinNO3anox ,
  kinO2anox  =  kinO2anox  ,DispO2    =  DispO2     ,
  DispNO3    =  DispNO3    ,DispNH3   =  DispNH3    ,
  DispODU    =  DispODU     )

# names of state variables and initial conditions
svarnames <- c("FDET", "SDET", "O2", "NO3", "NH3", "ODU")
nspec     <- length(svarnames)
Cini      <- rep(10, N*nspec)

#====================#
# budgetting         #
#====================#

IntegratedRate <- function(rate, depth = NULL)  {      # integrated rate for liquids
  if (is.null(depth))
    sum(rate* porGrid$mid * Grid$dx)
  else
    sum(rate* porGrid$mid * Grid$dx * (Grid$x.mid < depth))
}  
IntegratedRateSolid <- function(rate, depth = NULL) {  #                     solids
  if (is.null(depth))
    sum(rate* porGridSolid$mid * Grid$dx)
  else
    sum(rate* porGridSolid$mid * Grid$dx * (Grid$x.mid < depth))
}  

O2budget <- function(out) {

  flux <- out$O2flux - out$O2deepflux

  cons <- IntegratedRate(out$Nitri)*2 +
          IntegratedRate(out$OduOx)   +
          IntegratedRate(out$OxicMin)

  return(list(Flux = flux, Cons = cons, 
     Delta = flux - cons))
}

Nbudget <- function(out) {

  influx     <- out$NO3flux + out$NH3flux + out$Norgflux
  efflux     <- out$NO3deepflux + out$NH3deepflux + out$Norgdeepflux
  flux       <- influx - efflux
  
  NlossDenit <- IntegratedRate(out$Denitrific)*0.8
  Nadsorp    <- IntegratedRate(out$NH3adsorption)        

  return(list(Loss = NlossDenit, 
    influx = influx, burial = efflux, orgNflux = out$Norgflux,
    DINflux = out$NO3flux + out$NH3flux, 
    Denit = NlossDenit, Delta = flux - NlossDenit))
}

## =============================================================================
##
## Three runs of OMEXDIA, with different inputs of flux
##
## =============================================================================

pars1 <- pars            
pars1["MeanFlux"] <- 15000/12*100/365  # nmol/cm2/d - C deposition 

print(system.time(
 DIA1  <- steady.1D (y = Cini, func = OMEXDIA_C, names = svarnames,
                   parms = pars1, nspec = nspec, positive = TRUE)
))

pars2 <- pars
pars2["MeanFlux"] <- 25000/12*100/365  # nmol/cm2/d - C deposition 
print(system.time(
 DIA2  <- steady.1D (y = Cini, func = OMEXDIA_C, names = svarnames,
                   parms = pars2, nspec = nspec, positive = TRUE)
))

pars3 <- pars
pars3["MeanFlux"] <- 50000/12*100/365  # nmol/cm2/d - C deposition: 49gC/m2/yr
print(system.time(
 DIA3  <- steady.1D (y = DIA2$y, func = OMEXDIA_C, names = svarnames,
                   parms = pars3, nspec = nspec, positive = TRUE)
))

#====================#
# Plotting           #
#====================#

plot(DIA1, DIA2, DIA3,  
     which = c("O2", "NO3", "NH3", "ODU", "TOC"), 
     ylim = list(c(2, 0), c(5, 0), c(5,0), c(5,0), c(5,0)),
     grid = Grid$x.mid, lwd = 2, xlab = c(rep("mmol/m3",4), "%"),
     xyswap = TRUE, ylab = "depth, cm", 
     obspar = list(pch = ".", cex = 3))
     
plot.new()
CFlux <- c(pars1["MeanFlux"], pars2["MeanFlux"], pars3["MeanFlux"])*12/1e5*365
legtext <- paste(formatC(CFlux, 3), "gC/m2/yr")
legend ("bottom", col = 1:3, lty = 1:3, lwd = 2,
     legend = legtext, title = "import flux")

mtext(outer = TRUE, side = 3, line = -2, cex = 1.5, "OMEXDIA model")

#====================#
# budgetting         #
#====================#

Oxygen <- rbind(unlist(O2budget(DIA1)),  unlist(O2budget(DIA2)),
                unlist(O2budget(DIA3)))
print(Oxygen)
                

Nitrogen <- rbind(unlist(Nbudget(DIA1)), unlist(Nbudget(DIA2)),
                  unlist(Nbudget(DIA3)))
                 
print(Nitrogen)

