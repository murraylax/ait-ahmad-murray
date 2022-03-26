%**************************************************************************************************
%                      AIT Paper                                                  *
%                      Yamin Ahmad, UW-Whitewater                                                 *
%                      James Murray, UW-La Crosse                                                  *
%                      Initial Ver: 1.0 (Feb 2022)                                               *
%                      Current Ver: 1.0 (Feb 2022)                                               *
%**************************************************************************************************

% Variables are expressed in logs, except for parameters

var         
            x,              % Output gap
            infl,           % Current inflation
            r,              % Nominal int rate
            exg,            % Expected output gap
            exinfl,         % Expected future inflation
            pia,            % AIT 
            pib,            % Backwards looking inflation
            pif,            % Forwards looking inflation
            xix,            % Demand pull shocks
            xip;            % Cost push shocks
        
    
varexo      epsxt,          % Shock to IS Curve
            epspt,          % Shock to Phillips Curve
            epsmp;          % Shock to monetary policy rule



parameters  beta, sigma, gamma, kappa, lambda, rhox, rhop, rhor, psix, psip, deltab, deltaf, rstar, trials, per;
            load parameterfile1;

            beta = 0.99;        % Discount factor
            sigma = 2;          % Inverse Intertemporal elasticity
            kappa = 0.1;        % Marginal cost coefficient
            rhor = 0.7;         % Monetary policy persistence 
            psix = 0.5;         % CB weight on output gap 
            psip = 1.5;         % CB weight on inflation 
            rhox = 0.7;         % Persistence of demand shock
            rhop = 0.7;         % Persistence of inflation shock
            rstar = 0.02;       % Natural rate of interest
             
            % Key AIT Parameters loaded by parameterfile1
            
            % deltab = 0.25;      % Weight on Backwards looking inflation term 
            set_param_value('deltab',modparams(1));
            % deltaf = 0.25;      % Weight on Forwards looking inflation term
            set_param_value('deltaf',modparams(2));
            % gamma = 0;          % AIT weight on past inflation
            set_param_value('gamma',modparams(3));
            % lambda = 0;         % Weight on naive expectations/Minnesota Prior
            set_param_value('lambda',modparams(4));
            
            % Simulation parameters
            trials = 1000;      % Number of data-generating trials
            per = 1000;         % Number of periods


model;

    
    %*****************************    Definitional Equations and MARKET CLEARING    **************************
    
    x = exg - (r - exinfl - rstar)/sigma + xix;
    infl = beta * exinfl + kappa*x + xip;

    exg = lambda*x + (1-lambda)*x(+1);
    exinfl = lambda*infl + (1-lambda)*infl(+1);
    r = rhor*r(-1) + (1-rhor)*(psip*pia + psix*x) + epsmp;
    
    %**********************************    AIT equations    ***************************************

    pia = gamma*pib + (1-gamma)*pif;
    % pib = deltab*(...sum...((1-deltab)^(j))*infl(...-j...));
    pib = (1-deltab)*infl + (1-deltab)*pib(-1);
    % pif = deltaf*(...sum...((1-deltaf)^(j))*infl(...j...)); 
    pif = deltaf*infl(+1)+ (1-deltaf)*pif(+1);   % Weight on inflation has to be greater than or equal to 1 otherwise model is non-stationary and there is indeterminacy
    
    %*************************************    SHOCK Processes    ***************************************
    xix  = rhox * xix(-1) + epsxt;
    xip  = rhop * xip(-1) + epspt;
    % ximp  = ximp(-1) + epsmp;
       

end;

initval;

    xix = 0;
    xip = 0;
    % ximp = 0;
    r  = rstar;
    
    % NUMERICAL SOLUTIONS, GIVEN ABOVE PARAMETER VALUES
    x = 0;
    infl = 0.02;
    pib = 0.02;
    % pif = 0.02;
    
   % ANALYTICAL SOLUTIONS FOR THE REST OF THE VARIABLES
    pia = infl;
    exg = lambda*x;
    exinfl = lambda*infl;

   
end;

shocks;

var epsxt = 0.25;
var epspt = 0.25;
var epsmp = 0.01;

end;

% Determine steady state
% ==============================
% steady;

% Business cycle properties, theoretical moments
% ==============================================
check;
% model_info;
% model_diagnostics;

% Value of Shock put in directly in previous Dynare version.
% Sigma_e = [0.000049];
     %     epszt                         

% Run Stochastic Simulation
% ==============================
% stoch_simul(irf=20,periods = 5500, drop=500, hp_filter=1600,ORDER=1) x, infl, exg, exinfl, r;
