############################################################################
##
#F  Mean ( <L> )
##
##  <Description>
##
##     Returns the mean value for the list <A>L</A> of real numbers
##
##  </Description>
##
Mean := function(list)
   return Sum(list)/Size(list);
end;

############################################################################
##
#F  _Variance ( <L> )
##
##  <Description>
##
##     Returns the variance of the list <A>L</A> of real numbers
##
##  </Description>
##
_Variance := function(list)
  local Sxx, n;

  n := Size(list);
  Sxx := list*list-n*(Mean(list)^2);

  return Sxx/n; #some texts use n-1 instead of n
end;


############################################################################
##
#F  Covariance ( <xList>, <yList> )
##
##  <Description>
##
##     Returns the covariance for the list <A>xList</A> and <A>yList</A> of real numbers
##
##  </Description>
##
Covariance := function(xList, yList)
  local n, Sxy;

  n := Size(xList);
  if n <> Size(yList) then
    Error("the size of xList is not the same that yList ");
  fi;

  Sxy:= xList*yList-(n*Mean(xList)*Mean(yList));

  return Sxy/n;
end;

############################################################################
##
#F  LinearRegression ( <xList>, <yList> )
##
##  <Description>
##
##     Returns a record with properties a and b, that results from computing the
##     coefficients a and b of the linear regression
##
##         f(x)=a+bx,
##
##     using <A>xList</A> as the x_i values and <A>yList</A> as the f(x_i) values for the
##     linear regression.
##
##  </Description>
##
LinearRegression := function(xList, yList)
  local a,b, covXY, varX, varY, rSquared;

  covXY := Covariance(xList, yList);
  varX := Variance(xList);
  varY := Variance(yList);

  b := covXY/varX;
  a := Mean(yList)-b*Mean(xList);
  rSquared := (covXY^2)/(varX*varY);

  return rec(a:=a, b:=b, rSquared:=rSquared);
end;

############################################################################
##
#F  LogarithmicRegression ( <xList>, <yList> )
##
##  <Description>
##
##     Returns a record with properties a and b, that results from computing the
##     coefficients a and b of the logarithmic regression
##
##         f(x)=a \cdot b^x (which is equivalent to Ln(f(x))=Ln(a)+Ln(b)*x,
##
##     using <A>xList</A> as the x_i values and <A>yList</A> as the f(x_i) values for the
##     linear regression.
##
##  </Description>
##
LogarithmicRegression := function(xList, yList)
  local xListFloat, yListFloat, r;

  xListFloat := List(xList, Float);
  yListFloat := List(yList, Float);
  r := LinearRegression(xListFloat, List(yListFloat,Log));

  r.a := Exp(r.a);
  r.b := Exp(r.b);

  #Print("log regression:",r,"\n");

  return r;
end;
