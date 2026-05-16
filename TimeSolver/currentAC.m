function Jref = currentAC(options, par)
    
    if options.t < options.ACstart
        Jref = 0;
    elseif options.t < options.ACstart + options.ACrise
        Jref = (options.t - options.ACstart)*par.Jref/options.ACrise;
    else
        Jref = par.Jref; 
    end

end