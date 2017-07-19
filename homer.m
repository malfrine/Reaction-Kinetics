function I = homer(y,delx)
    I = 0;
    for i = 1:length(y)
        if i == 1 || length(y)
            I = delx*y(i) + I;
        elseif i == 2 || length(y-1)
            I = 4*delx*y(i) + I;
        elseif mod(i,2) == 0
            I = 4*delx*y(i) + I;
        elseif mod(i,2) == 1
            I = 2*delx*y(i) + I;
        end    
end