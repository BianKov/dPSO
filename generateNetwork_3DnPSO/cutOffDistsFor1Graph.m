function Rvalues = cutOffDistsFor1Graph(d,zeta,N,m,beta,T)


    Rvalues = zeros(N,1); %initialization
    for t=1:N
        if m<t-1 %if t-1<=m, then Rvalues(t) remains 0 (there is no node sampling)
            R_t_0 = 0; %initial value
            Rvalues(t) = fzero(@(x) ZERO_FCN(double(d),double(beta),double(zeta),double(T),double(t),x,double(m)), R_t_0);
        end
    end


end
