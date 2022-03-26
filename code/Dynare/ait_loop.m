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

% Set parameter values
deltafs = (0.01:0.01:1)';
deltabs = 0.25;
gammas = 0;
lambdas = 0;
neig=0;
resindet = zeros(length(deltafs),1);

% Loop over Dynare to check if parameter values generate indeterminacy
first_time = 1;
for i = 1:length(deltafs)
    modparams = [deltabs; deltafs(i); gammas; lambdas];

    % set_param_value('deltaf',deltafs(i));
    save parameterfile1 modparams;
    
    if first_time
        dynare ait_v1;
        E = oo_.dr.eigval;
        neig = numel(E(abs(E)>1));
        first_time=0;
    else
        set_param_value('deltaf',modparams(2,1));
        info = check(M_,options_,oo_);
        E = info;
        neig = numel(E(abs(E)>1));
    end
    

    if neig<3
        resindet(i)=1;  % Rank condition not satisfied - indeterminacy
    else
        resindet(i)=0;  % Rank condition satisfied
    end
end

clc;

% Report Results
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
