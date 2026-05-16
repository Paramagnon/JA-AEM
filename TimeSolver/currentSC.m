function Isc = currentSC(options,pol)
    
    if options.t < options.trise
        Isc = pol*options.Iset*options.t/options.trise;
    else
        Isc = pol*options.Iset; 
    end

end
