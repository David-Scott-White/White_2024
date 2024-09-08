function printDwellFit(dwellFits)
% David S. White 
% 2021-12-29 

% Mono exponential Fit 
tau0_value = dwellFits.monoExpTau; 
%tau0_error = dwellFits.monoExpTau - dwellFits.monoExpCI(1);
tau0_error = dwellFits.monoExpTau - dwellFits.monoExpCI(1);

% double exponetial --
tau1_value = dwellFits.biExpTau(1); 
%tau1_error = dwellFits.biExpTau(1) - dwellFits.biExpTauCI(1,1);
tau1_error = dwellFits.biExpTauSE(1);

tau2_value = dwellFits.biExpTau(2); 
%tau2_error = dwellFits.biExpTau(2) - dwellFits.biExpTauCI(1,2);
tau2_error = dwellFits.biExpTauSE(2);

amp1_value = dwellFits.biExpAmp(1); 
%amp1_error = dwellFits.biExpAmp(1) - dwellFits.biExpAmpCI(1,1);
amp1_error = dwellFits.biExpAmpSE(1);

amp2_value = dwellFits.biExpAmp(2); 
%amp2_error = dwellFits.biExpAmp(2) - dwellFits.biExpAmpCI(1,2);
amp2_error = dwellFits.biExpAmpSE(2);

disp('Dwell Summary: ')
disp(['>> N: ', num2str(dwellFits.nDwells)]); 
disp(['>> tau0: ', num2str(round(tau0_value,2)), ' ± ',  num2str(round(tau0_error,2))]);
disp(['>> tau1: ', num2str(round(tau1_value,2)), ' ± ',  num2str(round(tau1_error,2))]);
disp(['>> Amp1: ', num2str(round(amp1_value,2)), ' ± ',  num2str(round(amp1_error,2))]);
disp(['>> tau2: ', num2str(round(tau2_value,2)), ' ± ',  num2str(round(tau2_error,2))]);
disp(['>> Amp2: ', num2str(round(amp2_value,2)), ' ± ',  num2str(round(amp2_error,2))]);
disp(['>> LLR: ' , num2str(dwellFits.LLR(3))])
disp(['>> BIC: ' , num2str(dwellFits.BIC(3))])

end