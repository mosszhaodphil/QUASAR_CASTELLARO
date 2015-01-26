




function h=plot_ci(t,dm,p_low,p_median,p_high,color_bande,color_median)
ci = prctile(dm,[p_low p_median p_high],2);
    
    hf=fill([t' ; flipud(t')],[ci(:,1) ; flipud(ci(:,3))],color_bande);
    set(hf,'EdgeColor',color_median)
    hold on
    plot(t,ci(:,2),'LineWidth',2,'Color',color_median);


end




