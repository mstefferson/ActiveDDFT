%%
num_rec = size(C_rec, 3);
for ii = 1:num_rec
  cTemp = C_rec(:,:,ii);
  disp(min(cTemp))
end
%%
xind = 1;
yind = 1;
[~,~,nm] = size(Den_rec(:,:,:,1));
figure
ax = gca;
ax.YLim = [0 1];
for ii = 1:num_rec
  fTemp = reshape( Den_rec(xind, yind,:,ii), [1 nm] );
  plot(fTemp)
  disp( min(fTemp) )
  ax.YLim = [0 1];
  keyboard
end
%%
xind = 1;
phi_ind = 1;
[~,~,nm] = size(Den_rec(:,:,:,1));
figure
ax = gca;
ax.YLim = [0 1];
for ii = 1:num_rec
  fTemp = reshape( Den_rec(xind, :, phi_ind, ii), [1 nm] );
  plot(fTemp)
  disp( min(fTemp) )
  ax.YLim = [0 1];
  keyboard
end