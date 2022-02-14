function plotdelay(f0, s4p0, fignum)

    if (f0(1) == 0)
	f   = f0;
	s4p = s4p0;
    else
	if (f0(1) < 1e6)
	    df  = diff(f0(1:2));
	else
	    df = f0(1);
	end
	nf  = round(f0(end)/df);
	f   = df * [0:nf]';
	s4p = interpsnp(f0, s4p0, f);
    end

    s4pm = s4p2s4pm(s4p);

    [pdnl4p,  pd4p,  phnl4p,  ph4p ] = snp2pdnoloss(f, s4p);
    [pdnl4pm, pd4pm, phnl4pm, ph4pm] = snp2pdnoloss(f, s4pm);

    figure(fignum+1);
    hold off;
    plot(f/1e9, 20*log10(abs(s4pm(:,2,1))), 'o-', 'Disp', 'S21dd');
    hold on;
    plot(f/1e9, 20*log10(abs(s4pm(:,4,3))), 'o-', 'Disp', 'S21cc');
    plot(f/1e9, 20*log10(abs(s4p(:,2,1))), 'Disp', 'S21');
    plot(f/1e9, 20*log10(abs(s4p(:,4,3))), 'Disp', 'S43');
    plot(f/1e9, 20*log10(abs(s4p(:,4,1))), 'Disp', 'S41');
    plot(f/1e9, 20*log10(abs(s4p(:,2,3))), 'Disp', 'S23');
    set(gca, 'FontSize', 16);
    xlim([0 max(f)/1e9]);
    %ylim([5.0 11.0]);
    legend({},'FontSize',12);
    legend show;
    legend('Location','SouthWest');
    title('S21 Gain');
    xlabel('frequency (GHz)');
    ylabel('gain (dB)');
    grid on;

    figure(fignum+2);
    hold off;
    %plot(f(2:end)/1e9, - diff(ph4p(:,2,1)) ./ diff(f) / 2 / pi / 1e-9, 'Disp', 'GD21');
    %plot(f(2:end)/1e9, - diff(ph4p(:,4,3)) ./ diff(f) / 2 / pi / 1e-9, 'y', 'Disp', 'GD43');
    plot(f(2:end)/1e9, pd4p(2:end,2,1)   / 1e-9, 'Disp', 'PD21');
    hold on;
    plot(f(2:end)/1e9, pd4p(2:end,4,3)   / 1e-9, 'Disp', 'PD43');
    plot(f(2:end)/1e9, pdnl4p(2:end,2,1) / 1e-9, 'Disp', 'LCPD21');
    plot(f(2:end)/1e9, pdnl4p(2:end,4,3) / 1e-9, 'Disp', 'LCPD43');
    set(gca, 'FontSize', 16);
    xlim([0 max(f)/1e9]);
    %ylim([7.75 7.95]);
    legend({},'FontSize',12);
    legend show;
    legend('Location','SouthEast');
    title('Single-end Response Delay');
    xlabel('frequency (GHz)');
    ylabel('delay (ns)');
    grid on;

    figure(fignum+3);
    hold off;
    %plot(f(2:end)/1e9, - diff(ph4pm(:,2,1)) ./ diff(f) / 2 / pi / 1e-9, 'Disp', 'GD21dd');
    %plot(f(2:end)/1e9, - diff(ph4pm(:,4,3)) ./ diff(f) / 2 / pi / 1e-9, 'y', 'Disp', 'GD21cc');
    plot(f(2:end)/1e9, pd4pm(2:end,2,1)   / 1e-9, 'Disp', 'PD21dd');
    hold on;
    plot(f(2:end)/1e9, pd4pm(2:end,4,3)   / 1e-9, 'Disp', 'PD21cc');
    plot(f(2:end)/1e9, pdnl4pm(2:end,2,1) / 1e-9, 'Disp', 'LCPD21dd');
    plot(f(2:end)/1e9, pdnl4pm(2:end,4,3) / 1e-9, 'Disp', 'LCPD21cc');
    set(gca, 'FontSize', 16);
    xlim([0 max(f)/1e9]);
    %ylim([5.0 9.0]);
    %ylim([7.65 7.95]);
    legend({},'FontSize',12);
    legend show;
    legend('Location','SouthEast');
    title('Mixed-mode Response Delay');
    xlabel('frequency (GHz)');
    ylabel('delay (ns)');
    grid on;

    figure(fignum+5);
    hold off;
    plot(f(2:end)/1e9,((pdnl4p(2:end,2,1) - pdnl4p(2:end,4,3)))/1e-12, 'Disp', 'S21 intra-pair skew');
    hold on;
    plot(f(2:end)/1e9,((pdnl4p(2:end,1,2) - pdnl4p(2:end,3,4)))/1e-12, 'Disp', 'S12 intra-pair skew');
    set(gca, 'FontSize', 16);
    xlim([0 max(f)/1e9]);
    %ylim([-20 5]);
    legend({},'FontSize',12);
    legend show;
    legend('Location','SouthEast');
    title('Intra-pair skew');
    xlabel('frequency (GHz)');
    ylabel('LCPD21(12) - LCPD43(34) (ps)');
    grid on;

    figure(fignum+6);
    hold off;
    plot(f(2:end)/1e9,((pdnl4pm(2:end,2,1) - pdnl4pm(2:end,4,3)))/1e-12, 'Disp', 'S21 inter-mode skew');
    hold on;
    plot(f(2:end)/1e9,((pdnl4pm(2:end,1,2) - pdnl4pm(2:end,3,4)))/1e-12, 'Disp', 'S12 inter-mode skew');
    set(gca, 'FontSize', 16);
    xlim([0 max(f)/1e9]);
    %ylim([20 90]);
    legend show;
    legend('Location','SouthEast');
    title('Inter-mode skew');
    xlabel('frequency (GHz)');
    ylabel('LCPD21(12)dd - LCPD21(12)cc (ps)');
    grid on;

    drawnow();
    scr = get(0,'ScreenSize');
    pos = get(gcf, 'OuterPosition');
    figure(fignum+1); set(gcf, 'OuterPosition', [pos(3)*0 scr(4)-pos(4)*1 pos(3) pos(4)]);
    figure(fignum+2); set(gcf, 'OuterPosition', [pos(3)*1 scr(4)-pos(4)*1 pos(3) pos(4)]);
    figure(fignum+3); set(gcf, 'OuterPosition', [pos(3)*2 scr(4)-pos(4)*1 pos(3) pos(4)]);
    %figure(fignum+4); set(gcf, 'OuterPosition', [pos(3)*0 scr(4)-pos(4)*2 pos(3) pos(4)]);
    figure(fignum+5); set(gcf, 'OuterPosition', [pos(3)*1 scr(4)-pos(4)*2 pos(3) pos(4)]);
    figure(fignum+6); set(gcf, 'OuterPosition', [pos(3)*2 scr(4)-pos(4)*2 pos(3) pos(4)]);

end

