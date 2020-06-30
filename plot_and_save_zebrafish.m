function plot_and_save_zebrafish(domain_matrix,domain_matrix_ir, domain_matrix_X,xb_on,stage,type,listx)
close all
h=figure();
colorspec = {[0.5843, 0.8157, 0.9882],[242,242,242]/255, [255,204,0]/255, [255,239,251]/255};
% 



if strcmp(type,'ablation')==0
    [sizex_xi, sizey_xi]=size(domain_matrix_ir);
    [sizex_xm, sizey_xm]=size(domain_matrix);
    
    domain_matrix_ir=domain_matrix_ir(1:sizex_xi,min(sizey_xi,ceil(1*sizex_xi)):min(sizey_xi,ceil(2.2*sizex_xi)));
    domain_matrix_X=domain_matrix_X(1:sizex_xi,min(sizey_xi,ceil(1*sizex_xi)):min(sizey_xi,ceil(2.2*sizex_xi)));
    domain_matrix=domain_matrix(1:sizex_xm,min(sizey_xm,ceil(1*sizex_xm)):min(sizey_xi,ceil(2.2*sizex_xm)));
    
else
    domain_matrix=domain_matrix(1:end,51:110);
    domain_matrix_X=domain_matrix_X(1:end,102:220);
    domain_matrix_ir=domain_matrix_ir(1:end,102:220);
    
end


[sizex_xi, sizey_xi]=size(domain_matrix_ir); 
[sizex_xm, sizey_xm]=size(domain_matrix);
% 
% rectangle('Position',[0 0 sizey_xi  sizex_xi], 'FaceColor',colorspec{4},'EdgeColor','k',...
%     'LineWidth',6)
hold on

[y_id, x_id]=find(domain_matrix_ir==6);
[y_il, x_il]=find(domain_matrix_ir==5);
[y_xan, x_xan]=find(domain_matrix_X==4);
[y_xb, x_xb]=find(domain_matrix_X==2);

[y_mel, x_mel]=find(domain_matrix==1);

y_id=y_id+0.5*rand(length(y_id),1)-0.75;
x_id=x_id+0.5*rand(length(x_id),1)-0.75;
y_il=y_il+0.5*rand(length(y_il),1)-0.75;
x_il=x_il+0.5*rand(length(x_il),1)-0.75;
y_xan=y_xan+0.5*rand(length(y_xan),1)-0.75;
x_xan=x_xan+0.5*rand(length(x_xan),1)-0.75;
y_il=y_il+0.5*rand(length(y_il),1)-0.75;
x_il=x_il+0.5*rand(length(x_il),1)-0.75;
y_xb=y_xb+0.5*rand(length(y_xb),1)-0.75;
x_xb=x_xb+0.5*rand(length(x_xb),1)-0.75;

y_mel=2*y_mel;
x_mel=2*x_mel;
y_mel=y_mel+1.5*rand(length(y_mel),1)-1.25;
x_mel=x_mel+1.5*rand(length(x_mel),1)-1.25;
hold on



for i=1:length(y_il)
    if y_il(i)>sizex_xi/2
        plot(x_il(i),y_il(i),'^','MarkerFaceColor', colorspec{1},'MarkerEdgeColor', colorspec{1}, 'MarkerSize',1)
    else
        plot(x_il(i),y_il(i),'v','MarkerFaceColor',colorspec{1}, 'MarkerEdgeColor',colorspec{1},'MarkerSize',1)
    end
end
hold on

plot(x_id,y_id,'o','MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',2)
hold on
plot(x_xan,y_xan,'o','MarkerFaceColor',colorspec{3}, 'MarkerEdgeColor', colorspec{3}, 'MarkerSize',1) %[ 0.9100 0.4100 0.1700]
hold on
plot(x_mel,y_mel,'ko','MarkerFaceColor','k','MarkerSize',2)
hold on

if xb_on==1 
plot(x_xb,y_xb,'ro','MarkerFaceColor','r','MarkerSize',3)
end

% if 21+t<30
%     title([' \fontsize{30}PB ~ ',num2str(21+t),'dpf ']);
% elseif 21+t<39
%     title([' \fontsize{30}PR ~ ',num2str(21+t), 'dpf ']);
% elseif 21+t<44
%     title([' \fontsize{30}SP ~',num2str(21+t), 'dpf ']);
% elseif 21+t<51
%     title([' \fontsize{30}SR ~',num2str(21+t), 'dpf ']);
% elseif 21+t<71
%     title([' \fontsize{30}J ~',num2str(21+t), 'dpf ']);
% else
%     title([' \fontsize{30}J+ ~',num2str(21+t), 'dpf ']);
% end
% 
daspect([1 1 1])
xlim([0, sizey_xi]);
ylim([0, sizex_xi]);

hold on

if listx~=0
    plot([1,sizey_xi],[median(listx), median(listx)], 'b', 'LineWidth', 4)
end

t=text(1,sizex_xi+sizex_xi/10+5,stage);
t.FontSize=20;
t.FontWeight='bold';
% set(gca,'XTick',0:50:sizey_xi+2);
% set(gca,'YTick',0:50:sizex_xi+2);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'fontsize',20)
 rectangle('Position',[0 0 sizey_xi  sizex_xi], 'FaceColor','None','EdgeColor','k',...
     'LineWidth',1)
hold on

if strcmp(stage,'7.5mm') || strcmp(stage,'9mm') || strcmp(stage,'10.2mm') || strcmp(stage,'11.6mm')
    plot([3,26],[3,3], 'r', 'LineWidth', 3)
else
    plot([3,15],[3,3], 'r', 'LineWidth', 3)
end


if strcmp(stage,'7.5mm')
%     saveas(h,['Final_WT_fish_simulation/',type,'750','.fig'])
%     saveas(h,['Final_WT_fish_simulation/',type,'750','.png'])
    saveas(h,['Final_WT_fish_simulation/',type,'750','.eps'])
elseif strcmp(stage,'9mm')
%     saveas(h,['Final_WT_fish_simulation/',type,'900','.fig'])
%     saveas(h,['Final_WT_fish_simulation/',type,'900','.png'])
     saveas(h,['Final_WT_fish_simulation/',type,'900','.eps'])
elseif strcmp(stage,'10.2mm')
%     saveas(h,['Final_WT_fish_simulation/',type,'1020','.fig'])
%     saveas(h,['Final_WT_fish_simulation/',type,'1020','.png'])
    saveas(h,['Final_WT_fish_simulation/',type,'1020','.eps'])
elseif strcmp(stage,'11.6mm')
%         saveas(h,['Final_WT_fish_simulation/',type,'1160','.fig'])
%     saveas(h,['Final_WT_fish_simulation/',type,'1160','.png'])
    saveas(h,['Final_WT_fish_simulation/',type,'1160','.eps'])
else
    saveas(h,['simulations/',type,stage,'.eps'])
    saveas(h,['simulations/',type,stage,'.png'])
% saveas(h,['Final_WT_fish_simulation/',type,stage,'.fig'])
% saveas(h,['Final_WT_fish_simulation/',type,stage,'.png'])
end


