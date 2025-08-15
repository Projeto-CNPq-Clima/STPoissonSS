#' Vehicle theft events and neighborhood data for Natal/RN (Brazil)
#'
#' A dataset for spatio-temporal analysis of vehicle thefts in the city of Natal (RN), Brazil.
#' It contains three components:
#' 1. A matrix of event times (by neighborhood)
#' 2. UTM coordinates of neighborhood centroids
#' 3. Neighborhood covariates
#'
#' @section Time origin and units:
#' Event times were converted from strings of the form \code{"%d/%m/%Y %H:%M:%S"}
#' using \code{strptime} in R and expressed as numeric times referenced to
#' \code{"01/01/2012 00:00:00"}. The numeric values represent **seconds since**
#' \code{2012-01-01 00:00:00} (local time). To recover calendar datetimes:
#' \preformatted{
#'   origin <- as.POSIXct("2012-01-01 00:00:00", tz = "America/Fortaleza")
#'   dt <- origin + natal_theft_times[, 1]  # POSIXct for the first neighborhood
#' }
#'
#' @section Neighborhood order (columns/rows alignment):
#' The **columns of \code{natal_theft_times}**, the **rows of \code{natal_sites_utm}**
#' and the **rows of \code{natal_covariates}** follow exactly the order below.
#' Neighborhood names are kept exactly as provided.
#'
#' \strong{South Zone}
#' 1. Lagoa Nova; 2. Nova Descoberta; 3. Candelária; 4. Capim Macio; 5. Pitimbu;
#' 6. Neópolis; 7. Ponta Negra
#'
#' \strong{East Zone}
#' 8. Santos Reis; 9. Rocas; 10. Ribeira; 11. Praia do Meio; 12. Petrópolis;
#' 13. Areia Preta; 14. Cidade Alta; 15. Mãe Luiza; 16. Alecrim; 17. Barro Vermelho;
#' 18. Tirol; 19. Lagoa Seca
#'
#' \strong{West Zone}
#' 20. Cidade da Esperança; 21. Quintas; 22. Dix-Sept Rosado; 23. Bom Pastor;
#' 24. Nossa Senhora de Nazaré; 25. Felipe Camarão; 26. Cidade Nova; 27. Guarapes; 28. Planalto
#'
#' \strong{North Zone}
#' 29. Igapó; 30. Potengi; 31. Nossa Senhora da Apresentação; 32. Lagoa Azul;
#' 33. Pajuçara; 34. Redinha
#'
#' @format The \code{.RData} file contains three objects:
#' \describe{
#'   \item{\code{natal_theft_times}}{A numeric matrix of size \eqn{T \times 34}.
#'         Each column corresponds to one neighborhood (in the order above) and
#'         contains the event times (in seconds since \code{"2012-01-01 00:00:00"}).
#'         Rows represent individual events; columns may have different numbers of
#'         non-missing entries depending on neighborhood event counts.}
#'
#'   \item{\code{natal_sites_utm}}{A numeric matrix of size \eqn{34 \times 2}
#'         with UTM coordinates (Eastings, Northings) of the centroids of the 34
#'         neighborhoods, in the exact same order as above.
#'         Columns: \code{Easting}, \code{Northing}.}
#'
#'   \item{\code{natal_covariates}}{A numeric matrix of size \eqn{34 \times 3}
#'         with neighborhood-level covariates, one row per neighborhood in the same order.
#'         Columns (in order): \code{density} (population density),
#'         \code{pm_stations} (number of Polícia Militar stations),
#'         \code{squares} (number of public squares).}
#' }
#'
#' @details
#' Original raw event timestamps were provided as strings like
#' \code{"01/01/2012 00:00:00"} and parsed with
#' \preformatted{strptime(dados[i,], "%d/%m/%Y %H:%M:%S")}.
#' Neighborhood names are preserved verbatim (no translation).
#'
#' @note
#' The three objects were named in English for package consistency:
#' \itemize{
#'   \item Proposed file name: \code{natal_vehicle_theft_data.RData}
#'   \item Proposed object names: \code{natal_theft_times}, \code{natal_sites_utm}, \code{natal_covariates}
#' }
#' If you are migrating from legacy names, these correspond to:
#\code{Tempo_Ocorrencia_roubo} -> \code{natal_theft_times};
#\code{sites} -> \code{natal_sites_utm};
#\code{Covariáveis} -> \code{natal_covariates}.
#'
#' @examples
#' # Load the dataset included in the package
#' data(natal_vehicle_theft_data)
#'
#' # Check what is inside
#' names(natal_vehicle_theft_data)
#'
#' # Check the dimensions of each component
#' dim(natal_vehicle_theft_data$natal_theft_times)   # T x 34 (T = number of events)
#' dim(natal_vehicle_theft_data$natal_sites_utm)     # 34 x 2
#' dim(natal_vehicle_theft_data$natal_covariates)    # 34 x 3
#'
#' # Recover POSIXct datetimes for a given neighborhood (e.g., Lagoa Nova = column 1)
#' origin <- as.POSIXct("2012-01-01 00:00:00", tz = "America/Fortaleza")
#' theft_dt_lagoa_nova <- origin + natal_theft_times[, 1]
#'
#' # Example: build a simple sf object of centroids (if you know the UTM zone/CRS)
#' # library(sf)
#' # utm_crs <- 31985  # example only; set the correct SIRGAS 2000 / UTM zone for Natal
#' # pts <- st_as_sf(as.data.frame(natal_sites_utm), coords = c("Easting", "Northing"), crs = utm_crs)
#'
#' # Merge covariates with coordinates
#' # neigh_df <- cbind(neighborhood = c("Lagoa Nova","Nova Descoberta", ... "Redinha"),
#' #                   natal_covariates, natal_sites_utm)
#' }
#'
#' @keywords datasets spatial spatio-temporal crime point-process
#' @docType data
#' @name natal_vehicle_theft_data
#' @aliases natal_theft_times natal_sites_utm natal_covariates
#' @usage data(natal_vehicle_theft_data)
NULL
