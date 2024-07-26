x_sd = logspace(-3,3,10000);
y_sd = 0.0025./(x_sd+0.009);
loglog(x_sd, y_sd)

hold on;
x_sdm = logspace(-3,3,10000);
y_sdm = 0.004066./(x_sdm+0.01426);
loglog(x_sdm, y_sdm)

hold on;
x_sm = logspace(-3,3,10000);
y_sm = 0.00747./(x_sm+0.03531);
loglog(x_sm, y_sm)

hold on;
x_bcd = logspace(-3,3,10000);
y_bcd = 0.003176./(x_bcd+0.0061);
loglog(x_bcd, y_bcd)

%solar = char(9737)
xlabel('Radius (kpc)')

solar=char(9737)
ylabel(['Density (M_',solar, '/pc^{3})']);

%ylabel('Density (M_{\odot}/pc^{3}', 'Interpreter', 'latex')
title('Average SPARC density curves')


legend('Sd', 'Sdm', 'Sm', 'BCD')