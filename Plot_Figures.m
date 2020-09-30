if (PLOT_E)
    % Plot first few Even Energy levels
    plot(EvenIt,EvenE);
    title_txt = ['Energy vs Iterations (EVEN)'];
    title(title_txt);
    xlabel('Iterations[It]');
    ylabel('Energy[E]');
    
    % Plot first few Odd Energy levels
    figure;
    plot(OddIt,OddE);
    title_txt = ['Energy vs Iterations (ODD)'];
    title(title_txt);
    xlabel('Iterations[It]');
    ylabel('Energy[E]');
end

if (PLOT_T)
    % Plot Entropy
    figure;
    semilogx(temp(2:ITMAX),simp(2:ITMAX));
    title_txt = 'Entropy vs Temperature';
    title(title_txt);
    xlabel('Temperature[T]');
    ylabel('Entropy[Savg]');
    
    % Plot TChi
    figure;
    semilogx(temp(2:ITMAX),tchiimp(2:ITMAX));
    title_txt = 'TChi vs Temperature';
    title(title_txt);
    xlabel('Temperature[T]');
    ylabel('[TChi]');
    
    % Plot Specific Heat
    figure;
    semilogx(temp(2:ITMAX),Cvimp(2:ITMAX));
    title_txt = 'Specific Heat vs Temperature';
    title(title_txt);
    xlabel('Temperature[T]');
    ylabel('Specific Heat [Cv]');
end