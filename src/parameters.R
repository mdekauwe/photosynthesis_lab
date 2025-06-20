#' Parameters for the Farquhar Photosynthesis Model
#'
#' This script defines all the constants and physiological parameters used in 
#' the Farquhar photosynthesis model, including temperature dependencies, enzyme 
#' kinetics, and leaf optical properties. These are collected into a list (`p`) 
#' at the bottom of the script for convenient passing to model functions.
#'
#' Sources include Bernacchi et al. (2001), Medlyn et al. (2002), and others as 
#' noted.
#' 
#' @author Martin De Kauwe
#' @date 20/06/2025
#' @name parameters
NULL

Oi <- 210.0 # intercellular concentration of O2 [mmol mol-1]

# -------------------------------------------------------------------------
# Enzyme kinetic constants at 25 °C
# -------------------------------------------------------------------------
gamstar25 <- 42.75 # CO2 compensation point - base rate at 25 deg C / 298 K 
                   # [umol mol-1]
Kc25 <- 404.9      # Michaelis-Menten coefficents for carboxylation by Rubisco 
                   # at 25degC [umol mol-1] or 298 K
Ko25 <- 278.4      # Michaelis-Menten coefficents for oxygenation by Rubisco at
                   # 25degC [mmol mol-1]. Note value in Bernacchie 2001 is in 
                   #  mmol!! or 298 K
Vcmax25 <- 40.0          # max rate of rubisco activity at 25 deg or 298 K
Jmax25 <- Vcmax25 * 1.67 # potential rate of electron transport at 25 deg or 298 K


Rd25 <- 1.4 # Respiration rate at the reference temperature 25 deg C or 298 K 
# [deg K]
# Aspinwall et al 2016 New Phyt (WTC3 estimates)
#Rd25 = 1.4*0.7 # DK reduce Rday 30

Q10 <- 2.0 # ratio of respiration at a given temperature divided by respiration
           # at a temperature 10 degrees lower

# -------------------------------------------------------------------------
# Activation energies 
# -------------------------------------------------------------------------
Ec <- 79430.0   # for carboxylation [J mol-1]
Eo <- 36380.0   # for oxygenation [J mol-1]
Eag <- 37830.0  # at CO2 compensation point [J mol-1]
Eaj <-43790.0   # for Jmax [J mol-1] #30000.0 
Eav <- 51560.0  # for Vcmax [J mol-1] #60000.0 

# -------------------------------------------------------------------------
# Deactivation energies 
# -------------------------------------------------------------------------
Hdv <- 200000.0 # for Vcmax [J mol-1]
Hdj <- 200000.0 # for Jmax [J mol-1]

deltaSj <- 644.4338 #650.0 # entropy factor [J mol-1 K-1)
deltaSv <- 629.26 #650.0 # entropy factor [J mol-1 K-1)

theta_hyperbol <- 0.9995 # Smooth transition between two limitations
theta_J <- 0.7       # Curvature of the light response (-)
quantum_yield <- 0.3
tau <- 0.1           # leaf transmissivity [-] (VIS: 0.07 - 0.15)
                     # ENF: 0.05; EBF: 0.05; DBF: 0.05; C3G: 0.070
refl <- 0.1          # leaf reflectance [-] (VIS:0.07 - 0.15)
                     # ENF: 0.062;EBF: 0.076;DBF: 0.092; C3G: 0.11
absorptance <- 1.0 - refl - tau # unitless

# Leaf quantum yield (initial slope of the A-light response curve) [mol mol-1]
# this value (0.246) is slightly lower than MAESTRA (0.26)
# reflectance and transmittance should change to MAESTRA

#alpha = quantum_yield * absorptance # (Medlyn et al 2002)
alpha <- 0.24

# -------------------------------------------------------------------------
# Stomatal params
# -------------------------------------------------------------------------
g0 <- 1E-09 # residual stomatal conductance as net assimilation rate reaches 
            # zero (mol m-2 s-1), effectively zero
g1 <- 4.0 # slope of the sensitivity of stomatal conductance to assimilation
          # (mol m-2 s-1)

# -------------------------------------------------------------------------
# Leaf optical properties
# -------------------------------------------------------------------------
leaf_width <- 0.02 # leaf width (m)

# Cambell & Norman, 11.5, pg 178
# The solar absorptivities of leaves (-0.5) from Table 11.4 (Gates, 1980)
# with canopies (~0.8) from Table 11.2 reveals a surprising difference.
# The higher absorptivity of canopies arises because of multiple reflections
# among leaves in a canopy and depends on the architecture of the canopy.
SW_abs <- 0.8 # use canopy absorptance of solar radiation
emissivity_leaf <- 0.96 # leaf emissivity (-), Table 3, Wang and Leuning, 1998
emissivity_soil <- 0.94 # soil emissivity (-), Table 3, Wang and Leuning, 1998
soil_reflectance <- 0.1 #(same as MAESTRA) # Table 3, Wang and Leuning, 1998
k <- 0.5 # light extinction coefficient

# -------------------------------------------------------------------------
# Bundle key parameters into a list for use in model functions
# -------------------------------------------------------------------------
p <- list(Kc25=Kc25, Ko25=Ko25, Eo=Eo, Ec=Ec, Oi=Oi, gamstar25=gamstar25,
          Eag=Eag, Vcmax25=Vcmax25, Eav=Eav, deltaSv=deltaSv, Hdv=Hdv,
          Jmax25=Jmax25, Eaj=Eaj, deltaSj=deltaSj, Hdj=Hdj, theta_J=theta_J,
          alpha=alpha, g0=g0, g1=g1)
