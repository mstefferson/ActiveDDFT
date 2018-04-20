%% Plot of external potential test for thesis
saveMe = 1;
saveName = 'external_pot_validation';
%%
fileId = 'Hr_disks_quadV_1.00_2.00_diag1_N128321_ls105_bc0.016_fD00_ICgauss_p_SM6_t01_01';
dirloc = [ 'analyzedfiles/potValid/' fileId '/'];
filename = ['run_' fileId '.mat'];
load([dirloc filename])
filename = ['params_' fileId '.mat'];
load([dirloc filename])
%%
nt = timeObj.N_rec;
n1 = systemObj.n1;
n2 = systemObj.n2;
l2 = systemObj.l2;
l1 = systemObj.l1;
x1 = (l1/n1 * [-n1/2:n1/2-1]).';
x2 = (l2/n2 * [-n2/2:n2/2-1]).';
rho_slice = l2*reshape(Den_rec(:,1,1,:), [n1, nt]);
%% theory
k = particleObj.externalV{1}{2}(2);
t = linspace(0.1,2.1,9);
alpha_sig = 1 / k * (1 - exp(-2*k*t));
x = (10/n1 * [-n1/2:n1/2-1]).';
rho_theory = 1 ./ sqrt(2*pi*alpha_sig) .*...
  exp( -x.^2 ./ (2*alpha_sig) );
%%
figure()
fig = gcf;
fig.WindowStyle = 'normal';
fig.Position = [292 181 550 491];
hold on
ind = [1 2 4];
my_colors =  viridis(3);
my_leg = cell(2*length(ind),1);
for ii = 1:length(ind)
  id = ind(ii);
  p = area(x, rho_slice(:,id));
  p.FaceColor =  my_colors(ii,:);
  p.EdgeColor =  my_colors(ii,:);
  alpha(0.3)
  my_leg{2*ii-1} = ['DDFT $$t=$$' num2str(t(ii)-t(1))];
%   p.Color =  my_colors(1,:);
  p = plot(x, rho_theory(:,id),'--');
  p.Color =  my_colors(ii,:);
  my_leg{2*ii} = ['theory $$t=$$' num2str(t(ii)-t(1))];
%   p.Color =  my_colors(2,:);
end
ax = gca;
axis square
xlabel('Position $$x$$');
ylabel('Density profile $$\rho(x,t)$$');
hl = legend(my_leg);
hl.Interpreter = 'latex';
box on
hl.Position = [0.5811 0.6357 0.2850 0.2772];
if saveMe
  saveas(gcf, saveName,'png')
end