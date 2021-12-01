function IntegralResult_2 = Integral_2(d,beta,zeta,T,t,R_t)


    fun = @(x) Integral_1(d,beta,zeta,T,t,x,R_t);

    IntegralResult_2 = integral(fun,1,t);


end




