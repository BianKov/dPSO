function I = numIntegralInAnalyticCutoffDist(dim,Temp)

    integrand = @(q,d,T) (sin(q/2).^(-1./T)).*(sin(q).^(d-2));
    I = integral(@(q)integrand(q,double(dim),double(Temp)),0,pi,'AbsTol',1e-12,'RelTol',1e-25);

end
