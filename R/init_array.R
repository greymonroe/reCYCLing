#' Initialise a tandem repeat array
#'
#' Creates the starting array of \code{k0} identical repeat units used to seed
#' a simulation. The ancestral monomer can be drawn randomly from an alphabet
#' or supplied directly.
#'
#' @param l Integer. Length in base pairs of one monomer unit. Default 178.
#' @param k0 Integer. Number of initial copies (all identical). Default 10.
#' @param init_sequence_type Character. \code{"random"} (draw a new monomer
#'   from \code{alphabet}) or \code{"given"} (use \code{ancestor_seq}).
#' @param ancestor_seq Character scalar. Required when
#'   \code{init_sequence_type = "given"}; must have \code{nchar() == l}.
#' @param alphabet Character vector of allowed bases. Default
#'   \code{c("A","C","G","T")}.
#'
#' @return A character vector of length \code{k0}, every element the same
#'   monomer sequence.
#' @export
#'
#' @examples
#' # Random ancestral monomer, 10 copies of 30 bp
#' arr <- init_array(l = 30, k0 = 10)
#' length(arr)   # 10
#' nchar(arr[1]) # 30
#'
#' # Supply your own starting sequence
#' ancestor <- paste(rep("A", 30), collapse = "")
#' arr2 <- init_array(l = 30, k0 = 5,
#'                    init_sequence_type = "given",
#'                    ancestor_seq = ancestor)
init_array <- function(l      = 178,
                       k0     = 10,
                       init_sequence_type = c("random", "given"),
                       ancestor_seq       = NULL,
                       alphabet = c("A", "C", "G", "T")) {
  init_sequence_type <- match.arg(init_sequence_type)

  if (init_sequence_type == "given") {
    if (is.null(ancestor_seq))
      stop("ancestor_seq must be provided when init_sequence_type = 'given'")
    if (nchar(ancestor_seq) != l)
      stop("ancestor_seq must have length l = ", l)
    unit <- ancestor_seq
  } else {
    unit <- paste0(sample(alphabet, size = l, replace = TRUE), collapse = "")
  }

  rep(unit, k0)
}
