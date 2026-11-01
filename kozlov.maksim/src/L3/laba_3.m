clear; clc; close all;

K = 4.0;  T = 1.5;  T1 = 1.6;  T2 = 2.0;  T1a = 4.5;
s = tf('s');

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го', '5. Апериодическое 2-го', '6. Колебательное', ...
         '7. Консервативное', '8. Инерционное интегр.', '9. Реальное дифф.', ...
         '10. Форсирующее', '11. Изодромное'};

elements = @(K, T, T1, T2, T1a) {K/s, ...
                                 K*s, ...
                                 tf(K), ...
                                 K/(T*s + 1), ...
                                 K/(T2^2*s^2 + T1a*s + 1), ...
                                 K/(T2^2*s^2 + T1*s + 1), ...
                                 K/(T*s^2 + 1), ...
                                 K/(s*(T*s + 1)), ...
                                 K*s/(T*s + 1), ...
                                 K*(T*s + 1), ...
                                 K*(T*s + 1)/s};

W = elements(K, T, T1, T2, T1a);

w  = logspace(-2, 3, 2000);
wn = logspace(-1.3, 1.3, 4000);
w2 = logspace(-2, 2, 1000);

f1 = figure('Color', 'w', 'Position', [50 50 1400 900]);
f2 = figure('Color', 'w', 'Position', [80 80 1600 1000]);
f3 = figure('Color', 'w', 'Position', [110 110 1600 1000]);

for i = 1:numel(names)
    sys = W{i};

    [num, den] = tfdata(sys, 'v');
    Wa = @(x) polyval(num, 1i*x) ./ polyval(den, 1i*x);

    figure(f1); subplot(3,4,i); hold on; grid on;
    [re, im] = nyquist(sys, wn);
    re = squeeze(re).';  im = squeeze(im).';
    if isscalar(re), re = repmat(re, size(wn)); im = repmat(im, size(wn)); end
    Wv = Wa(wn);
    plot(re, im, 'b-', 'LineWidth', 2);
    plot(real(Wv), imag(Wv), 'r--', 'LineWidth', 1.5);
    plot(re(1),   im(1),   'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    plot(re(end), im(end), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    if max(re) - min(re) < 1e-3, xlim(round(mean(re)) + [-1 1]); end
    if i == 3, xlim([2 6]);   ylim([-1 1]); end
    if i == 7, xlim([-5 15]); ylim([-1 1]); end
    if i == 1, legend({'Числ.', 'Аналит.', '\omega=0.05', '\omega=20'}, 'Location', 'best'); end
    title(names{i}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Re'); ylabel('Im');

    [mag, ph] = bode(sys, w);
    Wv = Wa(w);
    figure(f2); subplot(3,4,i); hold on; grid on;
    semilogx(w, 20*log10(squeeze(mag)), 'b', 'LineWidth', 2);
    semilogx(w, 20*log10(abs(Wv)), 'r--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log');
    if i == 1, legend({'Числ.', 'Аналит.'}, 'Location', 'best'); end
    title(names{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('\omega, рад/с'); ylabel('L, дБ');

    phi_num = squeeze(ph).';
    phi_an  = rad2deg(angle(Wv));
    phi_an  = phi_an + 360*round((phi_num - phi_an)/360);
    figure(f3); subplot(3,4,i); hold on; grid on;
    semilogx(w, phi_num, 'b', 'LineWidth', 2);
    semilogx(w, phi_an, 'm--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log');
    if i == 1, legend({'Числ.', 'Аналит.'}, 'Location', 'best'); end
    title(names{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('\omega, рад/с'); ylabel('\psi, °');
end
exportgraphics(f1, '01_АФЧХ.png');
exportgraphics(f2, '02_ЛАЧХ.png');
exportgraphics(f3, '03_ЛФЧХ.png');

K_vals = [1 2 4 6 8];
K_lbl  = arrayfun(@(x) sprintf('K=%g', x), K_vals, 'UniformOutput', false);
WK = cell(1, numel(K_vals));
for j = 1:numel(K_vals)
    WK{j} = elements(K_vals(j), T, T1, T2, T1a);
end
f4 = plotSweep(WK, 1:numel(names), K_lbl, names, w2, 'L',   'Влияние K на ЛАЧХ', 3, 4);
exportgraphics(f4, '04_Влияние_K_ЛАЧХ.png');
f5 = plotSweep(WK, 1:numel(names), K_lbl, names, w2, 'phi', 'Влияние K на ЛФЧХ', 3, 4);
exportgraphics(f5, '05_Влияние_K_ЛФЧХ.png');

T_vals = [0.2 0.5 1 2 5];
T_lbl  = {'T/5', 'T/2', 'T', '2T', '5T'};
WT = cell(1, numel(T_vals));
for j = 1:numel(T_vals)
    m = T_vals(j);
    WT{j} = elements(K, T*m, T1*m, T2*m, T1a*m);
end
f6 = plotSweep(WT, 1:numel(names), T_lbl, names, w2, 'L',   'Влияние T на ЛАЧХ', 3, 4);
exportgraphics(f6, '06_Влияние_T_ЛАЧХ.png');
f7 = plotSweep(WT, 1:numel(names), T_lbl, names, w2, 'phi', 'Влияние T на ЛФЧХ', 3, 4);
exportgraphics(f7, '07_Влияние_T_ЛФЧХ.png');

function fig = plotSweep(Wsets, idx, labels, names, w, what, ttl, nr, nc)
    fig = figure('Color', 'w', 'Position', [30 30 1800 1000]);
    styles = {'-', '--', ':', '-.', '-'};
    widths = [3 2.5 2 2 1.5];

    for a = 1:numel(idx)
        i = idx(a);
        ax = subplot(nr, nc, a); hold on; grid on;

        for j = 1:numel(Wsets)
            [mag, ph] = bode(Wsets{j}{i}, w);
            if strcmp(what, 'L'), y = 20*log10(squeeze(mag)); else, y = squeeze(ph); end
            semilogx(ax, w, y, 'LineStyle', styles{j}, 'LineWidth', widths(j));
        end

        set(ax, 'XScale', 'log');
        title(ax, names{i}, 'FontSize', 10, 'FontWeight', 'bold');
        xlabel(ax, '\omega, рад/с');
        if strcmp(what, 'L')
            ylabel(ax, 'L, дБ');
        else
            ylabel(ax, '\psi, °');
            yl = ylim(ax);
            if diff(yl) < 20, c = mean(yl); ylim(ax, [c-90 c+90]); end
        end

        if a == 1, legend(ax, labels, 'Location', 'best'); end
    end
    sgtitle(ttl, 'FontWeight', 'bold');
end
