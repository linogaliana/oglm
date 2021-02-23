#' @export

print.oglmx <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...){

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  std <- diag(x$vcov)^0.5

  if (length(coef(x))){

    cat("Coefficients:\n")

    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co), co),
                                  1L, paste, collapse = "="), "]")

    formated_coef <- paste0(
      format(x$coefficients, digits = digits),
      signif_stars_vectorized(2*pnorm( -abs(x$coefficients/std)),
                              type = "none"),
      sprintf(" (%s)",format(std, digits = digits))
    )
    names(formated_coef) <- names(x$coefficients)

    print.default(
      formated_coef,
      print.gap = 2,
      quote = FALSE
    )

  } else{
    cat("No coefficients\n\n")
  }

  cat("\n")

  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")

  cat("Log likelihood:\t",
      format(signif(x$loglikelihood, digits)),
      "\nLog likelihood by observation:\t   ",
      format(signif(x$loglikelihood/attr(x$loglikelihood, "No.Obs"),
                    digits)),
      "\nAIC:\t", format(signif(AIC(x), digits)))
  cat("\n")
  invisible(x)

}

#' @export

print.oglmx.selection <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...){

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  std <- diag(x$vcov)^0.5
  std <- std[names(std) %in% names(x$estimate)]

  if (length(coef(x))){

    cat("Coefficients:\n")


    formated_coef <- paste0(
      format(x$estimate, digits = digits),
      signif_stars_vectorized(2*pnorm( -abs(x$estimate/std)),
                              type = "none"),
      sprintf(" (%s)",format(std, digits = digits))
    )
    names(formated_coef) <- names(x$estimate)

    print.default(
      formated_coef,
      print.gap = 2,
      quote = FALSE
    )

  } else{
    cat("No coefficients\n\n")
  }

  cat("\n")

  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")

  cat("Log likelihood:\t",
      format(signif(x$loglikelihood, digits)),
      "\nLog likelihood by observation:\t   ",
      format(signif(x$loglikelihood/attr(x$loglikelihood, "No.Obs"),
                    digits)),
      "\nAIC:\t", format(signif(AIC(x), digits)))
  cat("\n")
  invisible(x)

}
