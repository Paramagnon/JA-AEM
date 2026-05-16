function figureIEEE() 
    fig = figure;
    fig.Units               = 'centimeters';
    fig.Position(3)         = 8;
    fig.Position(4)         = 6;
    set(fig.Children, 'FontName', 'Times', 'FontSize', 9);
    set(gcf,'renderer','Painters')
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))    
    fig.PaperPositionMode   = 'auto';
end