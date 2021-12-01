function RES = ZERO_FCN(d,beta,zeta,T,t,R_t,m)


    RES = eta_(d).* Integral_2(d,beta,zeta,T,t,R_t)-m;


end




