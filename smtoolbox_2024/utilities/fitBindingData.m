%% Fit Binding Data
function [binding_fits,output] = fitBindingData(ligand_concentration, bound_probability, fit_options)

% binding_fits = zeros(length(fit_options), 3);
binding_fits = zeros(1, 3);
n = length(bound_probability);

% initial guesses for all
bmax0 = 1;
h0 = 1;
% kd0 = 1e-6; % (M)

[a,b] = min(abs(0.5-bound_probability));
kd0 =ligand_concentration(b);

disp(newline)
% Results returned as: [bmax, h, kd; ...]
switch fit_options
    case 1
        % give bmax, return kd, h
        hill = @(p, x, b) (b.*x.^p(1))./(x.^p(1)+p(2));
        costFunc = @(p) (sum((hill(p, ligand_concentration, bmax0) -  bound_probability).^2) / n).^0.5;
        p = fminsearch(costFunc, [h0,  kd0]);
        binding_fits(1) = bmax0;
        binding_fits(2) = p(1);
        binding_fits(3) = p(2);
        disp('Free Parm: h, kd');
        fprintf ('Bmax is %1.2f h is %1.2f kd is %1.2e\n',  binding_fits(:))
        
    case 2
        % give h, return kd bmax
        hill = @(p, x, h) (p(1).*x.^h)./(x.^h+p(2));
        costFunc = @(p) (sum((hill(p, ligand_concentration, h0) -  bound_probability).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0,  kd0]);
        binding_fits(1) = p(1);
        binding_fits(2) = h0;
        binding_fits(3) = p(2);
        disp('Free Parm: bmax, kd');
        fprintf ('Bmax is %1.2f h is %1.2f kd is %1.2e\n',  binding_fits(:))
        
    case 3
        % give bmax and h, return kd
        hill = @(p, x, b, h) (b.*x.^h)./(x.^h+p(1));
        costFunc = @(p) (sum((hill(p, ligand_concentration, bmax0, h0) -  bound_probability).^2) / n).^0.5;
        p = fminsearch(costFunc, kd0);
        binding_fits(1) = bmax0;
        binding_fits(2) = h0;
        binding_fits(3) = p(1);
        disp('Free Parm: kd');
        fprintf ('Bmax is %1.2f h is %1.2f kd is %1.2e\n',  binding_fits(:))
        
    case 4
        % Give nothing, return bmax, kd, h
        hill = @(p, x) (p(1).*x.^p(2))./(x.^p(2)+p(3));
        costFunc = @(p) (sum((hill(p, ligand_concentration) -  bound_probability).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, h0, kd0]);
        binding_fits(1) = p(1);
        binding_fits(2) = p(2);
        binding_fits(3) = p(3);
        disp('Free Parm: bmax, h, kd');
        fprintf ('Bmax is %1.2f h is %1.2f kd is %1.2e\n',  binding_fits(:))
        
    case 5
        % Give bmax and kd, return h 
        hill = @(p, x, b, kd) (b.*x.^p(1))./(x.^p(1)+kd);
        costFunc = @(p) (sum((hill(p, ligand_concentration, bmax0, kd0) -  bound_probability).^2) / n).^0.5;
        p = fminsearch(costFunc, h0);
        binding_fits(1) = bmax0;
        binding_fits(2) = p(1);
        binding_fits(3) = kd0;
        disp('Free Parm: h');
        fprintf ('Bmax is %1.2f h is %1.2f kd is %1.2e\n',  binding_fits(:))
        
end

%n = floor(log(abs(binding_fits(3)))./log(10));
%ligand_range = logspace(n-4, n+4,1000);
n1 = floor(log(abs(min(ligand_concentration)/10))./log(10));
n2 = floor(log(abs(max(ligand_concentration)*10))./log(10));
ligand_range = logspace(n1, n2,1000);

yhat = (binding_fits(1).*ligand_range.^binding_fits(2))./(ligand_range.^binding_fits(2)+binding_fits(3));
output = [ligand_range(:), yhat(:)];

end