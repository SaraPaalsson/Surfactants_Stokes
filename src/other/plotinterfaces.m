function plotinterfaces(z,idx,rho)
hold off
for j = 1:size(idx,1)
%     plot([z(idx(j,1):idx(j,2));z(idx(j,1))],'k')
    scatter(real([z(idx(j,1):idx(j,2));z(idx(j,1))]),imag([z(idx(j,1):idx(j,2));z(idx(j,1))]),...
        15,[rho(idx(j,1):idx(j,2));rho(idx(j,1))],'filled')
    hold on
    plot(z(idx(j,1):idx(j,2)),'k-')
    colorbar
end