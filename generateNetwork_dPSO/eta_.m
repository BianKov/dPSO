function eta = eta_(d)


    if (mod(d,2) == 0) && (d~=0) %d is even       
        eta = (factorial(d/2-1).*factorial((d-2)/2).*(2.^(d-2))) ./ (factorial(d-2).*pi);
    elseif (mod(d,2) == 1) && (d~=1) %d is odd
        eta = (factorial(d-1)) ./ (factorial((d-1)/2-1).*factorial((d-1)/2).*(2.^(d-1)));
    else %error
        disp('----- ERROR -----');
        disp(['d = ' num2str(d) ' is not possible']);
        eta = NaN;
    end    


end
