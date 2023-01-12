function plotgrid(figure_data,x,y)
hold on
    %along y axis
    if y == 1
        for i = [0.5:1:size(figure_data,2)]
            plot( [i i],get(gca,'ylim'),'k','linestyle','-','LineWidth',0.0002);
        end
    end
    
    if x == 1
        %along x axis
        for i= [0.5:1:size(figure_data,1)]
            plot(get(gca,'xlim'), [i i],'k','linestyle','-','LineWidth',0.0002);
        end
    end
end