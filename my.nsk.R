my.nsk <<- function(x,kn,bk,tf,lag=NULL){
  if(tf=="integral" | tf=="avg"){
    out <- splines2::nsk(x,knots=kn,Boundary.knots=bk,integral=TRUE)
  }
  if(tf=="current"){
    out <- splines2::nsk(x,knots=kn,Boundary.knots=bk)
  }
  if(tf == "cut_integral"){
    stopifnot(!is.null(lag))
    out <- splines2::nsk(x,knots=kn,Boundary.knots=Bk,integral=TRUE) - splines2::nsk(pmax(0,x-lag),knots=kn,Boundary.knots=bk,integral=TRUE)
  }
  if(tf=="deriv"){
    out <- splines2::nsk(x,knots=kn,Boundary.knots=Bk, derivs=1)
  }
  if(tf=="shared_only"){
    stop("shared-ranefs-only not implemented")
  }
  out
}
