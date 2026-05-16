function Jc = returnJcHTS(Top,B)
    % Source: https://gitlab.cern.ch/steam/steam-material-library/-/blob/master/models/c/general/CFUN_HTS_JcFit_SUPERPOWER_v1.c
    % Fitted constants - SUPERPOWER
    g0       = 0.039445688;     % [-]
    g1       = 0.276598888;     % [-]
    g2       = 0.047699303;     % [-]
    g3       = 0.055878401;     % [-]
    Tc0      = 93;              % [K]
    p_c      = 0.400002709;     % [-]
    q_c      = 2.940671164;     % [-]
    Bi0_c    = 140;             % [T]
    gamma_c  = 2.999950877;     % [-]
    alpha_c  = 2.99996E+12;     % [A*T/m^2]
    nu       = 1.720055527;     % [-]
    n0       = 1.169541488;     % [-]
    n1       = 1.599990183;     % [-]
    n2       = 1.131992777;     % [-]
    p_ab     = 0.825413689;     % [-]
    q_ab     = 5.495925721;     % [-]
    Bi0_ab   = 250;             % [T]
    a        = 2.05E-05;        % [-]
    gamma_ab = 4.998961387;     % [-]
    alpha_ab = 7.10415E+13;     % [A*T/m^2]
    thetaOff = -2.179824999;    % [deg]
    
    Bmin = 0.01; % [T]   
    Bnorm = max(Bmin,B);
        
    thetaIn = 0; % Effect of magnetic field angle not considered
    theta = deg2rad(abs(thetaIn - thetaOff));
    
    if (mod(theta,pi) <= pi/2)
        thetaMod = mod(theta,pi);
    else
        thetaMod = pi - mod(theta,pi);
    end

    % Temperature scaling factors
    t_T  = Top/Tc0;
    t_n0 = 1-t_T^n0;
    t_n1 = 1-t_T^n1;

    % Irreversibility field  
	Bi_ab = Bi0_ab*(t_n1^n2+a*t_n0);
	Bi_c  = Bi0_c*t_n0;
    
	% Magnetic flux density	
	b_ab  = Bnorm / Bi_ab;
	b_c   = Bnorm / Bi_c;
	
    % Flags	
	TCrit_Flag   = t_T  >= 1 || t_T  < 0;
	babCrit_Flag = b_ab >= 1 || b_ab < 0;		
	bcCrit_Flag  = b_c  >= 1 || b_c  < 0;
	
    if (~TCrit_Flag)
        if (~babCrit_Flag)
            % critical current density ab-plane
            Jc_ab =	(alpha_ab/Bnorm*b_ab^p_ab) * (1-b_ab)^q_ab * (t_n1^n2 + a*t_n0)^gamma_ab;
        else
            % field over Bi
            Jc_ab = 0;
        end
        if (~bcCrit_Flag)
            % Critical current density c-plane
            Jc_c = (1/Bnorm)*alpha_c*(b_c^p_c)*((1-b_c)^q_c)*((1-(t_T^n0))^gamma_c);
        else
            % Field over Bi
            Jc_c = 0;
        end

        % Anisotropy factor
        g = g0 + g1*exp(-g2*Bnorm*exp(g3*Top));
        Jc = min(Jc_ab,Jc_c)+ max(0,(Jc_ab-Jc_c))/(1+((pi/2-thetaMod)/g)^nu);                
    else
        % temperature over Tc0
        Jc = 0;
    end
end