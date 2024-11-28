
functions {
  vector shapescale_to_rates(real shape, real scale,
			     int nprates,
			     vector train_data_x, // ntrain
			     matrix train_data_y, // ntrain x nprates 
			     matrix train_data_m, // ntrain x nprates
			     int spline
			     ){
    vector[nprates] canpars;
    for (i in 1:nprates){
      if (spline==1)
	canpars[i] = spline_interp_linear(shape, train_data_x, train_data_y[,i], train_data_m[,i]);
      else 
	canpars[i] = spline_interp_hermite(shape, train_data_x, train_data_y[,i], train_data_m[,i]);
    }
    vector[nprates] rates = canpars_to_rates(canpars, nprates) / scale;

    return rates;
  }

  vector canpars_to_rates(vector canpars,
			  int nprates){
    vector[nprates] ret;
    int nphase = (nprates+1) %/% 2;
    real qsoj1 = canpars[1];
    vector[nphase-1] incqsoj = canpars[2:nphase];
    vector[nphase-1] pabs_notlast = canpars[(nphase+1):(2*nphase-1)];
    vector[nphase]   qsoj = append_row(qsoj1, qsoj1 + cumulative_sum(incqsoj));
    vector[nphase-1] qsoj_notlast = qsoj[1:(nphase-1)];
    vector[nphase-1] arate_notlast = qsoj_notlast .* pabs_notlast;
    real arate_last = qsoj[nphase];
    vector[nphase]   arate = append_row(arate_notlast, arate_last);
    vector[nphase-1] prate = qsoj_notlast - arate_notlast;
    ret = append_row(prate, arate);
    return ret;
  }
  
  real spline_interp_linear(real x, // assumes inside x0
			    vector x0, // assumes increasing order
			    vector y0, vector m){
    real ret;
    int i = findinterval(x, x0);
    real dx = x0[i+1] - x0[i];
    real dy = y0[i+1] - y0[i];
    ret =  y0[i] + (x - x0[i]) * dy / dx; 
    return ret; 
  }

  real spline_interp_hermite(real x,
			     vector x0, // assumes increasing order
			     vector y0, vector m){
    real ret;
    int i = findinterval(x, x0);
    real h = x0[i+1] - x0[i];
    real t = (x - x0[i])/h; 
    real t1 = t - 1;
    real h01 = t*t*(3 - 2*t);
    real h00 = 1 - h01;
    real tt1 = t*t1;
    real h10 = tt1*t1;
    real h11 = tt1*t;
    ret = y0[i]*h00 + h*m[i]*h10 + y0[i+1]*h01 + h*m[i+1]*h11;
    return ret; 

  }

  // Assumes x is inside x0, will break if outside. 
  int findinterval(real x, vector x0){
    int i = 1;
    while ((i <= rows(x0)) && (x > x0[i+1])) {
      i = i+1;
    }
    return i;
  }

  
}
