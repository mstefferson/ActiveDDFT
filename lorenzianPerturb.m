function rho = lorenzianPerturb(rho, systemObj,perturbObj,gridObj)
% get perturbation
perturb = lorenzianPerturb(rho, systemObj, perturbObj, gridObj);
rho = rho + perturb;
end
