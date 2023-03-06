function plotresults(hfig,poles,results,marker)
%%PLOTRESULTS plot results of the matrix-function times vector 
% approximation   
figure(hfig)
% subplot(1,2,1)
% hold on
% semilogy(results.numpoles,results.abserrors,...
%     "DisplayName",poles.name,...
%     "LineWidth",2);
% hold off
% ylabel('Absolute Error')
% xlabel("Number of poles")
% subplot(1,2,2)
hold on
semilogy(results.numpoles,results.relerrors,...
    "DisplayName",poles.name,...
    "LineWidth",2,"Marker",marker);
hold off
ylabel('Relative Error')
xlabel("Number of poles")
legend('Location','northeast')
% subplot(1,2,1)
% set(gca,"YScale","log")
% axis square tight
% grid on
% subplot(1,2,2)
set(gca,"YScale","log")
axis square tight
grid on

end