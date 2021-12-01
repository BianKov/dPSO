function IntegralResult_1 = Integral_1(d,beta,zeta,T,t,s,R_t)

    if T<(1/(d-1))
        fun = @(x) (sin(x).^(d-2)) ./ (1+(sin(x/2).^(1./T)).*(s.^(beta./T)).*(t.^((2-beta)./T)).*exp(-zeta.*R_t./(2*T)));
    elseif (1/(d-1))<T
        fun = @(x) (sin(x).^(d-2)) ./ (1+(sin(x/2).^(1./T)).*(s.^(beta.*(d-1))).*(t.^((2-beta).*(d-1))).*exp(-zeta.*R_t./(2*T)));
    else
        'T can not be 1/(d-1)'
    end
        
    IntegralResult_1 = integral(fun,0,pi,'ArrayValued',true);


end




