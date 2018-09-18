function F = PLOTGPSHIP(XR,YR,YS,s2,names,figure_c)

% Plot results from simulation

% set(gcf,'Color','w');

[m,n]=size(YS);

for i=1:n
    FileName=['figure',num2str(figure_c)];
    f1=figure(figure_c);%figure('Name', names{i});
    set(0,'DefaultFigureWindowStyle','docked')
    plotgp(f1,XR,YR(:,i),YS(:,i), sqrt(s2(:,i)),names{i});
    if figure_c==1 ||figure_c==5||figure_c==9
        p=mtit('a)',...
	     'fontsize',22,...
	     'xoff',0,'yoff',0,'FontName','MS Reference Sans Serif');
    elseif figure_c==2 ||figure_c==6||figure_c==10
        p=mtit('b)',...
	     'fontsize',22,...
	     'xoff',-0.0,'yoff',0,'FontName','MS Reference Sans Serif');
    elseif figure_c==3 ||figure_c==7||figure_c==11
        p=mtit('c)',...
	     'fontsize',22,...
	     'xoff',0,'yoff',0,'FontName', 'MS Reference Sans Serif');
     elseif figure_c==4 ||figure_c==8||figure_c==12
        p=mtit('d)',...
	     'fontsize',22,...
	     'xoff',0,'yoff',0,'FontName', 'MS Reference Sans Serif');
     end
    
%     set(handle,'Position',[-100,0,0],'FontSize',18);
    saveas(f1,FileName,'meta');
    savefig(FileName);
    figure_c=figure_c+1;
end
F=figure_c;
end