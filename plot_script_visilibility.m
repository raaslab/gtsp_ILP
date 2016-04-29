clf; set(gcf,'position',[200 500 700 600]); hold on;
axis equal; axis on; axis([X_MIN X_MAX Y_MIN Y_MAX]);

patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
       'w' , 'linewidth' , 1.5 );
for i = 2 : size(environment,2)
    patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
           'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
end



plot3(guards(1:6,1), guards(1:6,2), 0.4*ones(1,length(guards(1:6,2))), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
text(guards_x,guards_y, 0.4*ones(1,length(guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);


plot3(targets_x, targets_y, 0.4*ones(1,length(targets_x)), 's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
text(targets_x,targets_y, 0.4*ones(1,length(targets_x)), target_nam,'HorizontalAlignment','left','FontSize',13);
