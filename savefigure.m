function savefigure(name,number,printfigs)
% Save figure(number) as f-name-number.pdf/.png if printfigs==true
%
ffname = sprintf('..//Figures//fig-%s-%i', name, number);
h=figure(number);
set(h,'Units','centimeters');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
if(printfigs==true)
    %savefig(sprintf('%s.fig',ffname));
    print(sprintf('%s.pdf',ffname),'-dpdf');
    %print(sprintf('%s.png',ffname),'-dpng');
end
end