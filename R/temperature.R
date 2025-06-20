
#' Arrhenius temperature response function
#'
#' @param k25 Rate parameter value at 25°C (298.15 K)
#' @param Ea Activation energy [J mol^-1]
#' @param Tk Leaf temperature [K]
#' @return Temperature-adjusted rate parameter
#' @export
arrh <- function(k25, Ea, Tk) {
  return ( k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk)) )
}

#' Peaked Arrhenius temperature response function
#'
#' @param k25 Rate parameter value at 25°C
#' @param Ea Activation energy [J mol^-1]
#' @param Tk Leaf temperature [K]
#' @param deltaS Entropy factor [J mol^-1 K^-1]
#' @param Hd Deactivation energy [J mol^-1]
#' @return Temperature-adjusted rate parameter
#' @export
peaked_arrh <- function(k25, Ea, Tk, deltaS, Hd) {
  arg1 <- arrh(k25, Ea, Tk)
  arg2 <- 1.0 + exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
  arg3 <- 1.0 + exp((Tk * deltaS - Hd) / (Tk * RGAS))
  
  return ( arg1 * arg2 / arg3 )
}

#' Michaelis-Menten constant for CO2/O2 at given temperature
#'
#' @param p Parameter list containing Kc25, Ko25, Ec, Eo, Oi
#' @param Tleaf Leaf temperature [K]
#' @return Michaelis-Menten constant Km [µmol mol^-1]
#' @export
calc_michaelis_menten_constants <- function(p, Tleaf) {
  
  Kc <- arrh(p$Kc25, p$Ec, Tleaf) # coeff for carboxylation, Kc (umol mol−1)
  Ko <- arrh(p$Ko25, p$Eo, Tleaf) # coeff for oxygenation, Ko (mmol mol−1)
  Km <- Kc * (1.0 + p$Oi / Ko)    # the Michaelis–Menten constant (umol mol−1)
  
  return ( Km )
}