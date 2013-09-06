% Plot the result for Kerker+Pulay, Pulay and Elliptic+Pulay

load KPE_Res_set1

figure(1)
set(0,'DefaultFigurePaperPosition', [0 0 6 6])
semilogy(VP{6},'r-');
hold on
semilogy(VKP{6}, 'b-');
semilogy(VEP{6}, 'k');
legend('Pulay','Kerker+Pulay','This work');
axis tight
xlabel('iter');
ylabel('log $(||r_k||/||V_k||)$');
