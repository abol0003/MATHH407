
load("workspacesteps.mat")


figure
h1 = plot(time(:), Data(:,3));
yticks(0:0.01:0.2)
ylabel("r(m)");
xlabel("Time (s)");


dt1 = datatip(h1, 9.98, Data(999,3), "Location", "southeast", "FontSize", 15);
dtend1 = datatip(h1, 12.2, Data(1221,3), "Location", "northeast", "FontSize", 15);

dt2 = datatip(h1, 19.98, Data(1999,3), "Location", "southeast", "FontSize", 15);
dtend2 = datatip(h1, 22.33, Data(2234,3), "Location", "northeast", "FontSize", 15);

dt3 = datatip(h1, 29.98, Data(2999,3), "Location", "southeast", "FontSize", 15);
dtend3 = datatip(h1, 32.34, Data(3235,3), "Location", "southeast", "FontSize", 15);

hold on

h2 = plot(time(:), r_desired(:));
legend("r (m)", "r_d (m)", "Location","northwest")
title("Step Responses", "FontSize", 18)

dtr0 = datatip(h2, 5, r_desired(500), "Location", "northwest", "FontSize", 15);
dtr1 = datatip(h2, 10, r_desired(1000), "Location", "northwest", "FontSize", 15);
dtr2 = datatip(h2, 20, r_desired(2000), "Location", "northwest", "FontSize", 15);
dtr3 = datatip(h2, 30, r_desired(3000), "Location", "northwest", "FontSize", 15);
