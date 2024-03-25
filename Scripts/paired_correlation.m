function [fig, ax,t]=paired_correlation(data,x_vars,y_vars, sig_level,savname,varargin)
% sig_level should be a scalar between 0 and 1 representing the desired significance level for color coding
% paired_correlation(tdata,col,col, 0.05,savname,...
%    20,'MarkerFaceColor',LC(2,:),'MarkerEdgeColor','none')

m=length(y_vars);
n=length(x_vars);

% Calculate aspect ratio
aspect_ratio = n / m;

% Get display size
screen_size = get(groot,'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);

% Calculate figure size based on aspect ratio and display size
fig_width = 0.5 * screen_width; % example value
fig_height = fig_width / aspect_ratio/1.3 ;

% Set figure position and size
fig_xpos = 0.1 * screen_width; % example value
fig_ypos = 0.1 * screen_height; % example value
fig_position = [fig_xpos, fig_ypos]; % calculated value
fig_size = [fig_width, fig_height]; % calculated value

fig = figure;
fig.Position = [fig_position, fig_size];
set(gcf,'Color','w')

t= tiledlayout(m,n,TileSpacing="compact",Padding="compact");
ax=[];
for i =1:m
    tax=[];
    for j =1:n
        tax=[tax,nexttile];
    end
    ax =[ax;tax];
end

for i =1:m
    for j =1:n
        x=data.(x_vars{j});y= data.(y_vars{i});

        if strcmpi(x_vars{j},y_vars{i})
            histogram(ax(i,j),x,Normalization="pdf",EdgeColor='none')
        else
            scatter(ax(i,j),x,y,varargin{:});hold(ax(i,j),'on')
            
            lm=fitlm(x,y);
            if lm.Coefficients.pValue(2) < sig_level
                color = 'r';
            else
                color = 'k';
            end
            xx = min(x):0.01:max(x);
            yy = predict(lm,xx');
            plot(ax(i,j),xx, yy, '-','LineWidth',1.5,Color= [1 1 1]*0.1)

            
            lh=text(ax(i,j),0.05,0.95,"\it R^2="+num2str(lm.Rsquared.Ordinary,2), 'Units','normalized');
            lh.Color=color;
            lh.FontSize=10;
%             ylim(ax(i,j),[0,inf]);xlim(ax(i,j),[0,inf]);
            tx=ax(i,j);
            ab = [0, min(tx.XLim(2), tx.YLim(2))];
            plot(tx,ab, ab, '--','Linewidth',0.5, Color=[1 1 1]*0.25)
        end
    end
end

for i = 1:m
    ylabel(ax(i,1),y_vars(i))
end

for  j =1:n
    xlabel(ax(m,j),x_vars(j))
end

if m>=1 && n>=2
    set(ax(1:m-1,:),'XTickLabel',[])
    set(ax(:,2:n),'YTickLabel',[])
end


set(ax, 'Fontsize', 9, 'LineWidth', 0.45,'box','on', ...
    'Xcolor', [1, 1, 1]*0.25, 'Ycolor', [1, 1, 1]*0.25,...
    'Xgrid','off','Ygrid','off')

% export_fig(fig,"results/"+savname,'-r300')

if ~isempty(savname)
    exportgraphics(fig,savname,Resolution=300)
end


end
