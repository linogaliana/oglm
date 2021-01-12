testthat::context("print")

set.seed(242)
n<-250
x1<-sample(c(0,1),n,replace=TRUE,prob=c(0.75,0.25))
x2<-vector("numeric",n)
x2[x1==0]<-sample(c(0,1),n-sum(x1==1),replace=TRUE,prob=c(2/3,1/3))
z<-rnorm(n,0.5)
# create latent outcome variable
latenty<-0.5+1.5*x1-0.5*x2+0.5*z+rnorm(n,sd=exp(0.5*x1-0.5*x2))
# observed y has four possible values: -1,0,1,2
# threshold values are: -0.5, 0.5, 1.5.
y<-vector("numeric",n)
y[latenty< -0.5]<--1
y[latenty>= -0.5 & latenty<0.5]<- 0
y[latenty>= 0.5 & latenty<1.5]<- 1
y[latenty>= 1.5]<- 2
dataset<-data.frame(y,x1,x2)


mod1 <- oglm::oglmx(y ~ x1 + x2 + z,
                    data=dataset,link="probit",constantMEAN=FALSE,
                    constantSD=FALSE,delta=0,threshparam=NULL)

output <- capture.output(
  print(mod1)
)

# print oglmx objects -------------

blanks <- which(output == "")
digits = max(3L, getOption("digits") - 3L)

# CALL PART ===============

testthat::expect_equal(
  capture.output(
    cat("\nCall:  ", paste(deparse(mod1$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  ),
  c("",as.character(output[(blanks[1]+1):(blanks[1]+2)]),"")
)


# COEFFICIENT PART ===========
std <- diag(mod1$vcov)^0.5


print_coef <- function(x){

  cat("Coefficients:\n")

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

}



testthat::expect_equal(
  capture.output(
    print_coef(mod1)
  ),
  c(as.character(output[(blanks[2]+1):(blanks[3]-1)]))
)


# PERFORMANCE PART ===============

print_perf <- function(x){

  cat("Log likelihood:\t",
      format(signif(x$loglikelihood, digits)),
      "\nLog likelihood by observation:\t   ",
      format(signif(x$loglikelihood/attr(x$loglikelihood, "No.Obs"),
                    digits)),
      "\nAIC:\t", format(signif(AIC(x), digits)))


}

testthat::expect_equal(
  capture.output(
    print_perf(mod1)
  ),
  c(as.character(output[(blanks[3]+1):length(output)]))
)
