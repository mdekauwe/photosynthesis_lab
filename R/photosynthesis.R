#' Model to calculate net leaf-level C3 photosynthesis following the Farquhar, 
#' von Caemmerer & Berry (1980) model of C3 photosynthesis.
#'
#' Photosynthesis is as the minimum of two limitations: Rubisco-limited
#' photosynthesis and RuBP-regeneration-limited photosynthesis. The model
#' assumes that Rubisco-limited photosynthesis follows a Michaelis–Menten
#' response function modified to account for a competitive inhibitor, oxygen.
#
#' This code is ignoring the third TPU or export limitation
#'
#'
#' @author Martin De Kauwe
#' @email mdekauwe@gmail.com
#' @date 20/06/2025
#' @name photosynthesis
NULL

#' Calculate net leaf-level C3 photosynthesis rate following Farquhar et al.
#' (1980)
#'
#' @param p List of model parameters
#' @param Tleaf Leaf temperature [K]
#' @param PAR Photosynthetically active radiation [µmol m^-2 s^-1]
#' @param Cs CO2 concentration at leaf surface [µmol mol^-1]
#' @param vpd Vapor pressure deficit at leaf surface [kPa]
#' @param peaked_Vcmax Logical, use peaked Arrhenius for Vcmax (default TRUE)
#' @param peaked_Jmax Logical, use peaked Arrhenius for Jmax (default TRUE)
#'
#' @return A list containing:
#' \item{An}{Net assimilation rate [µmol m^-2 s^-1]}
#' \item{Ac}{Rubisco-limited assimilation [µmol m^-2 s^-1]}
#' \item{Aj}{RuBP regeneration-limited assimilation [µmol m^-2 s^-1]}
#' \item{gsc}{Stomatal conductance to CO2 [mol m^-2 s^-1]}
#' \item{Vcmax}{Temperature adjusted Vcmax [µmol m^-2 s^-1]}
#' \item{Cic}{Intercellular CO2 when Rubisco-limited [µmol mol^-1]}
#' \item{Rd}{Day respiration rate [µmol m^-2 s^-1]}
#'
#' @export
calc_photosynthesis <-function(p, Tleaf, PAR, Cs, vpd, peaked_Vcmax=TRUE,
                               peaked_Jmax=TRUE) {

  #   References:
  #   ----------
  #   * Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical
  #     model of photosynthetic CO2 assimilation in leaves of C3 species.
  #     Planta, 149, 78-90.
  #   * Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley, P.C.,
  #     Kirschbaum, M.U.F., Leroux, X., Montpied, P., Strassemeyer, J.,
  #     Walcroft, A., Wang, K. and Loustau, D. (2002) Temperature response of
  #     Parameters of a biochemically based model of photosynthesis. II.
  #     A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
  #

  # Recycle vectors to same length
  max_len <- max(length(PAR), length(Cs), 1)
  PAR <- rep_len(PAR, max_len)
  Cs <- rep_len(Cs, max_len)
  Tleaf <- rep_len(Tleaf, max_len)
  
  # calculate temp dependancies of Michaelis-Menten constants for CO2, O2
  Km <- calc_michaelis_menten_constants(p, Tleaf)

  # CO2 compensation point in the absence of mitochondrial respiration
  # [umol mol-1]
  gamma_star <- arrh(p$gamstar25, p$Eag, Tleaf)

  # Calculate the maximum rate of Rubisco activity (Vcmax), accounting for
  # temperature dependancies
  if (peaked_Vcmax) {
    Vcmax <- peaked_arrh(p$Vcmax25, p$Eav, Tleaf, p$deltaSv, p$Hdv)
  } else {
    Vcmax <- arrh(p$Vcmax25, p$Eav, Tleaf)
  }

  # Calculate the potential rate of electron transport (Jmax), accounting for
  # temperature dependancies
  if (peaked_Jmax) {
    Jmax <- peaked_arrh(p$Jmax25, p$Eaj, Tleaf, p$deltaSj, p$Hdj)
  } else {
    Jmax <- arrh(p$Jmax25, p$Eaj, Tleaf)
  }

  # Leaf mitochondrial respiration in the light or day respiration
  # (umol m-2 s-1). Following Collatz et al. (1991), assume Rd 1.5% of Vcmax
  Rd <- 0.015 * Vcmax

  # Rate of electron transport, which is a function of absorbed PAR
  J <- calc_electron_transport_rate(p, PAR, Jmax)
  Vj <- J / 4.0

  gs_over_a <- calc_stomatal_coeff(p, Cs, vpd)

  # Solution when Rubisco activity is limiting
  Cic <- solve_ci(p, gs_over_a, Rd, Cs, gamma_star, Vcmax, Km)

  # Solution when electron transport rate is limiting
  Cij <- solve_ci(p, gs_over_a, Rd, Cs, gamma_star, Vj, 2.0*gamma_star)

  # Catch for low PAR, issue with Vj
  Cic <- ifelse(is_close(PAR, 0.0) | is_close(Vj, 0.0), Cs, Cic)
  Cij <- ifelse(is_close(PAR, 0.0) | is_close(Vj, 0.0), Cs, Cij)

  # Rate of photosynthesis when Rubisco activity is limiting
  Ac <- assim(Cic, gamma_star, Vcmax, Km)

  # Rate of photosynthesis when RuBP regeneration is limiting
  Aj <- assim(Cij, gamma_star, Vj, 2.0*gamma_star)

  # Catch for negative Ci and instances where Ci > Cs
  Ac <- ifelse(Cic <= 0.0 | Cic > Cs, 0.0, Ac)
  Aj <- ifelse(Cic <= 0.0 | Cic > Cs, 0.0, Aj)

  # When below light-compensation points, assume Ci=Ca.
  Aj <- ifelse(Aj <= Rd + 1E-09,
               assim(Cs, gamma_star, Vj, 2.0*gamma_star), Aj)

  # Hyperbolic minimum of Ac and Aj to smooth over discontinuity when moving
  # from electron # transport limited to rubisco limited photosynthesis
  A <- -mapply(quadratic, 1.0-1E-04, Ac+Aj, Ac*Aj, large=TRUE)

  # Net photosynthesis rate (umol m-2 s-1)
  An <- A - Rd

  # Calculate conductance to CO2 (mol m-2 s-1)
  gsc <- max(p$g0, p$g0 + gs_over_a * An)

  # Calculate conductance to water (mol m-2 s-1)
  gsw <- gsc * GSC_TO_GSW

  return ( list(An=An, Ac=Ac, Aj=Aj, gsc=gsc, Vcmax=Vcmax, Cic=Cic, Rd=Rd) )
}

#' Calculate assimilation rate given parameters
#'
#' @param Ci Intercellular CO2 concentration [µmol mol^-1]
#' @param gamma_star CO2 compensation point [µmol mol^-1]
#' @param a1 Rate parameter
#' @param a2 Rate parameter
#' @return Assimilation rate [µmol m^-2 s^-1]
#' @export
assim <- function(Ci, gamma_star, a1, a2) {
  return ( a1 * (Ci - gamma_star) / (a2 + Ci) )
}

#' Calculate electron transport rate for given absorbed PAR and Jmax
#'
#' @param p Parameter list containing theta_J and alpha
#' @param PAR Photosynthetically active radiation [µmol m^-2 s^-1]
#' @param Jmax Maximum electron transport rate [µmol m^-2 s^-1]
#' @return Electron transport rate J [µmol m^-2 s^-1]
#' @export
calc_electron_transport_rate <-function(p, PAR, Jmax) {
  A <- p$theta_J
  B <- -(p$alpha * PAR + Jmax)
  C <- p$alpha * PAR * Jmax
  J <- mapply(quadratic, A, B, C, large=FALSE)

  return (J)
}

#' Solve intercellular CO2 concentration Ci from quadratic (Leuning 1990)
#'
#' @param p Parameter list containing g0
#' @param gs_over_a Stomatal coefficient
#' @param Rd Day respiration [µmol m^-2 s^-1]
#' @param Cs Leaf surface CO2 concentration [µmol mol^-1]
#' @param gamma_star CO2 compensation point [µmol mol^-1]
#' @param gamma Parameter (Vcmax or Vj)
#' @param beta Parameter (Km or 2*gamma_star)
#' @return Intercellular CO2 concentration Ci [µmol mol^-1]
#' @export
solve_ci <- function(p, gs_over_a, Rd, Cs, gamma_star, gamma, beta) {
  A = p$g0 + gs_over_a * (gamma - Rd)

  arg1 <- (1. - Cs * gs_over_a) * (gamma - Rd)
  arg2 <- p$g0 * (beta - Cs)
  arg3 <- gs_over_a * (gamma * gamma_star + beta * Rd)
  B <- arg1 + arg2 - arg3

  arg1 <- -(1.0 - Cs * gs_over_a)
  arg2 <- (gamma * gamma_star + beta * Rd)
  arg3 <- p$g0 * beta * Cs
  C <- arg1 * arg2 - arg3

  Ci <- mapply(quadratic, A, B, C, large=TRUE)

  return ( Ci )
}

#' Calculate stomatal coefficient gs_over_a for Medlyn stomatal model
#'
#' @param p Parameter list containing g1, g0
#' @param Cs CO2 concentration at leaf surface [µmol mol^-1]
#' @param vpd Vapor pressure deficit [kPa]
#' @return Stomatal coefficient gs_over_a [mol m^-2 s^-1]
#' @export
calc_stomatal_coeff <- function(p, Cs, vpd) {

  vpd <- pmax(vpd, 0.05) # avoid very low VPD causing problems

  # 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
  # because we are calculating conductance to CO2!
  if (any((is_close(Cs, 0.0)))) {
    gs_over_a = 0.0
  } else {
    gs_over_a <- (1.0 + p$g1 / sqrt(vpd)) / Cs
  }

  return ( gs_over_a )
}
