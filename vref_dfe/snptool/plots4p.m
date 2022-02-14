function plots4p(f, s4p, fignum)

    s4pm = s4p2s4pm(s4p);

    %ui = 1/53.125e9;
    ui = 1/56e9;
    fmax = 1 / ( ui / 32 * 2);
    df = 10e6;
    nf = round(fmax/df);
    fi = [0:nf]'*df;

    s4pi = interpsnp(f, s4p, fi);

    s4pim = s4p2s4pm(s4pi);
    [i4pim,t4pim,d4pim] = snp2impl(fi, s4pim);
    u4pim = cumsum(i4pim);
    b21dd = u4pim(:,2,1) - interp1(t4pim(:,2,1), u4pim(:,2,1), t4pim(:,2,1) - ui, 'spline', u4pim(1,2,1));
    b21cc = u4pim(:,4,3) - interp1(t4pim(:,4,3), u4pim(:,4,3), t4pim(:,4,3) - ui, 'spline', u4pim(1,4,3));

    figure(fignum+1);
    plot(f/1e9, db([s4p(:,2,1) s4p(:,4,3) s4p(:,4,1) s4p(:,2,3) s4p(:,1,1) s4p(:,3,3) s4p(:,3,1) s4p(:,1,3)]), 'LineWidth', 2);
    set(gca,'FontSize',16);
    legend('S21', 'S43', 'S41', 'S23', 'S11', 'S33', 'S31', 'S13','Location','SouthWest');
    xlim([0 max(f)/1e9]);
    xlabel('Frequency (GHz)');
    ylabel('dB');
    title('Single-end Frequency Response');
    grid on;

    figure(fignum+2);
    %plot(f/1e9, db([s4pm(:,2,1) s4pm(:,4,1) s4pm(:,1,1) s4pm(:,2,2)]), 'LineWidth', 2);
    plot(f/1e9, db([s4pm(:,2,1) s4pm(:,1,1) s4pm(:,2,2)]), 'LineWidth', 2);
    set(gca,'FontSize',16);
    %legend('S21dd', 'S21cd', 'S11dd', 'S22dd', 'Location', 'SouthWest');
    legend('S21dd', 'S11dd', 'S22dd', 'Location', 'SouthWest');
    xlim([0 max(f)/1e9]);
    ylim([-60 0]);
    xlabel('Frequency (GHz)');
    ylabel('dB');
    title('Mixed-mode Frequency Response');
    grid on;

    figure(fignum+3);
    hold off;
    plot(t4pim(:,2,1)/1e-12, i4pim(:,2,1), '-', 'LineWidth', 2);
    hold on;
    plot(t4pim(:,4,3)/1e-12, i4pim(:,4,3), '-', 'LineWidth', 1);
    set(gca,'FontSize',16);
    legend('I21dd', 'I21cc');
    xlabel('Time (ps)');
    title('Thru Impulse Response');
    grid on;

    figure(fignum+4);
    hold off;
    plot(t4pim(:,2,1)/1e-12, u4pim(:,2,1), '-', 'LineWidth', 2);
    hold on;
    plot(t4pim(:,4,3)/1e-12, u4pim(:,4,3), '-', 'LineWidth', 1);
    set(gca,'FontSize',16);
    legend('U21dd','U21cc','Location','SouthEast');
    xlabel('Time (ps)');
    title('Thru Step Response');
    grid on;

    figure(fignum+5);
    hold off;
    plot(t4pim(:,1,1)/1e-12, u4pim(:,1,1), 'LineWidth', 2);
    hold on;
    plot(t4pim(:,2,2)/1e-12, u4pim(:,2,2), 'LineWidth', 1);
    set(gca,'FontSize',16);
    legend('U11dd', 'U22dd','Location','SouthEast');
    xlabel('Time (ps)');
    title('Reflection Step Response');
    grid on;

    figure(fignum+6);
    hold off;
    plot(t4pim(:,2,1)/1e-12, b21dd, '-', 'LineWidth', 2);
    hold on;
    plot(t4pim(:,4,3)/1e-12, b21cc, '-', 'LineWidth', 1);
    set(gca,'FontSize',16);
    legend('B21dd', 'B21cc');
    xlabel('Time (ps)');
    title('Thru Single-Bit Response');
    grid on;

    drawnow();
    scr = get(0,'ScreenSize');
    pos = get(gcf, 'OuterPosition');
    figure(fignum+1); set(gcf, 'OuterPosition', [pos(3)*0 scr(4)-pos(4)*1 pos(3) pos(4)]);
    figure(fignum+2); set(gcf, 'OuterPosition', [pos(3)*1 scr(4)-pos(4)*1 pos(3) pos(4)]);
    figure(fignum+3); set(gcf, 'OuterPosition', [pos(3)*2 scr(4)-pos(4)*1 pos(3) pos(4)]);
    figure(fignum+4); set(gcf, 'OuterPosition', [pos(3)*0 scr(4)-pos(4)*2 pos(3) pos(4)]);
    figure(fignum+5); set(gcf, 'OuterPosition', [pos(3)*1 scr(4)-pos(4)*2 pos(3) pos(4)]);
    figure(fignum+6); set(gcf, 'OuterPosition', [pos(3)*2 scr(4)-pos(4)*2 pos(3) pos(4)]);

end
