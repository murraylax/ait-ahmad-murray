% Program:          Report Results for Average Inflation Targeting Paper
%**************************************************************************************************
%                               AIT Paper                                                         *
%                               Yamin Ahmad, UW-Whitewater                                        *
%                               James Murray, UW-LaCrosse                                         *
%                               March 2022                                                        *
%**************************************************************************************************
%
% Program:  Find which values of Delta_f generate indeterminacy
% Date:     March 2022
% Current Version:  1.0 (March 2022) - Premiliminary Working Version
% Original Version: 1.0 (March 2022)- Preliminary working version


fprintf(1, 'Please choose from the following: \n');
fprintf(1, '0: Loop ONLY over Delta_f, otherwise Loop over Delta_f and: \n');
fprintf(1, '1: Loop over Delta_b \n');
fprintf(1, '2: Loop over Gamma \n');
fprintf(1, '3: Loop over Lambda \n');
modver = input('Please enter choice from above: ');
fprintf(1, 'Please enter a vector consisting of [lowerbound upperbound] for your choice above: \n');
fprintf(1, 'Enter values between (0 1] for your choice above');
bb = input('[Lowerbound Upperbound] =');
bb= bb(:);
lb = bb(1); ub = bb(2);
ngrid = 101;   % Default number of gridpoints

% Set run flags
pltgrph = 1;    % 1 - Plot graphs; 0 - do not plot (default)

% Set default parameter values
deltafs = (0.01:0.01:1)';
deltabs = 0.25;
gammas = 0;
lambdas = 0;

if modver==1
    deltabs = linspace(lb,ub,ngrid)';
    lvar = deltabs; lbl = '\delta_b';
elseif modver==2
    gammas = linspace(lb,ub,ngrid)';
    lvar = gammas;  lbl = '\gamma';
elseif modver==3
    lambdas = linspace(lb,ub,ngrid)';
    lvar = lambdas; lbl = '\lambda';
else
    deltafs = linspace(lb,ub,ngrid)';
    lvar = 1;
end

neig=0;
resindet = zeros(length(deltafs),length(deltabs));

tstart = tic;

% Loop over Dynare to check if parameter values generate indeterminacy
first_time = 1;
for j = 1:length(lvar)
    for i = 1:length(deltafs)

        if modver==1
            modparams = [lvar(j); deltafs(i); gammas; lambdas];
        elseif modver==2
            modparams = [deltabs; deltafs(i); lvar(j); lambdas];
        elseif modver==3
            modparams = [deltabs; deltafs(i); gammas; lvar(j)];
        else
            modparams = [deltabs; deltafs(i); gammas; lambdas];
        end

        % set_param_value('deltaf',deltafs(i));
        save parameterfile1 modparams;

        if first_time
            dynare ait_v1;
            E = oo_.dr.eigval;
            neig = numel(E(abs(E)>1));
            first_time=0;
        else
            set_param_value('deltaf',modparams(2,1));
            if modver==1
                set_param_value('deltab',modparams(1,1));
            elseif modver==2
                set_param_value('gamma',modparams(3,1));
            elseif modver==3
                set_param_value('lambda',modparams(4,1));
            end
            info = check(M_,options_,oo_);
            E = info;
            neig = numel(E(abs(E)>1));
        end


        if neig<3
            resindet(i,j)=1;  % Rank condition not satisfied - indeterminacy
        else
            resindet(i,j)=0;  % Rank condition satisfied
        end
    end
end

clc;

% Report Results

if length(lvar)==1  % For specific value of delta_b
    fprintf(1,'============================================ \n');
    fprintf(1,'Results \n');
    fprintf(1,'============================================ \n');
    fprintf(1,'Model parameterization: \n');
    fprintf(1,'Delta_b = %1.3f \n',deltabs);
    fprintf(1,'gamma = %1.3f \n',gammas);
    fprintf(1,'lambda = %1.3f \n',lambdas);
    varnames = {'Delta_f','Indeterminacy Flag'};
    restab = table(deltafs,resindet,'VariableNames',varnames);
    disp(restab);
    fprintf(1,'-------------------------------------------- \n');
    if pltgrph ==1
        figure; plot(deltafs,resindet); xlabel('\delta_f'); ylabel('Indeterminacy'); title (['Indeterminacy for Specific value of \delta_b =', num2str(deltabs)]);
    end

else
    fprintf(1,'============================================ \n');
    fprintf(1,'Results \n');
    fprintf(1,'============================================ \n');

    if pltgrph==1
        figure; 
        surf(deltafs,lvar,resindet');
        xlabel('\delta_f'); ylabel(lbl);
        title('Indterminacy Region');
    end
end

tstop = toc(tstart);
elaptime(tstop);