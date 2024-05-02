function plot_traces(time, GNSS_pos, GNSS_vel, DR_pos, DR_vel, ...
    DR_GNSS_pos, DR_GNSS_vel, heading)
Define_Constants;

%%  Position solution (2-by-2)
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(2, 2, 1);
plot(GNSS_pos(:, 2), GNSS_pos(:, 1), 'Color', '#4DBEEE');
ylabel('Latitude (°)');
xlabel('Longitude (°)');
title('GNSS');

subplot(2, 2, 2);
plot(DR_pos(:, 2) * rad_to_deg, DR_pos(:, 1) * rad_to_deg, 'Color', '#A2142F');
ylabel('Latitude (°)');
xlabel('Longitude (°)');
title('Dead Reckoning');

subplot(2, 2, 3);
plot(DR_GNSS_pos(:, 2) * rad_to_deg, DR_GNSS_pos(:, 1) * rad_to_deg, 'Color', '#EDB120');
ylabel('Latitude (°)');
xlabel('Longitude (°)');
title('DR/GNSS Integrated Navigation');

subplot(2, 2, 4);
p1 = plot(GNSS_pos(:, 2), GNSS_pos(:, 1), 'Color', '#4DBEEE');
hold on
p2 = plot(DR_pos(:, 2) * rad_to_deg, DR_pos(:, 1) * rad_to_deg, 'Color', '#A2142F');
p3 = plot(DR_GNSS_pos(:, 2) * rad_to_deg, ...
    DR_GNSS_pos(:, 1) * rad_to_deg, 'Color', '#EDB120');
ylabel('Latitude (°)');
xlabel('Longitude (°)');
lgd = legend([p1, p2, p3], {'GNSS', 'Dead Reckoning', 'DR/GNSS Integrated Navigation'});
lgd.FontSize = 8;
title('All');
hold off
saveas(gcf, 'output/position.eps', 'epsc');



%% Velocity graph (2-by-2)
% (north)
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(2, 2, 1);
plot(time, GNSS_vel(:, 1), 'Color', '#4DBEEE');
xlabel('Time (s)');
ylabel('Velocity north (m/s)');
title('GNSS');

subplot(2, 2, 2);
plot(time, DR_vel(:, 1), 'Color', '#A2142F');
xlabel('Time (s)');
ylabel('Velocity north (m/s)');
title('Dead Reckoning');

subplot(2, 2, 3);
plot(time, DR_GNSS_vel(:, 1), 'Color', '#EDB120');
xlabel('Time (s)');
ylabel('Velocity north (m/s)');
title('DR/GNSS Integrated Navigation');

subplot(2, 2, 4);
p1 = plot(time, GNSS_vel(:, 1), 'Color', '#4DBEEE');
hold on
p2 = plot(time, DR_vel(:, 1), 'Color', '#A2142F');
p3 = plot(time, DR_GNSS_vel(:, 1), 'Color', '#EDB120');
xlabel('Time (s)');
ylabel('Velocity north (m/s)');
lgd = legend([p1, p2, p3], {'GNSS', 'Dead Reckoning', 'DR/GNSS Integrated Navigation'});
lgd.FontSize = 8;
hold off
title('All');
saveas(gcf, 'output/velocity_north.eps', 'epsc');

%% Velocity graph (2-by-2)
% (east)
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(2, 2, 1);
plot(time, GNSS_vel(:, 2), 'Color', '#4DBEEE');
xlabel('Time (s)');
ylabel('Velocity east (m/s)');
title('GNSS');

subplot(2, 2, 2);
plot(time, DR_vel(:, 2), 'Color', '#A2142F');
xlabel('Time (s)');
ylabel('Velocity east (m/s)');
title('Dead Reckoning');

subplot(2, 2, 3);
plot(time, DR_GNSS_vel(:, 2), 'Color', '#EDB120');
xlabel('Time (s)');
ylabel('Velocity east (m/s)');
title('DR/GNSS Integrated Navigation');

subplot(2, 2, 4);
p1 = plot(time, GNSS_vel(:, 2), 'Color', '#4DBEEE');
hold on
p2 = plot(time, DR_vel(:, 2), 'Color', '#A2142F');
p3 = plot(time, DR_GNSS_vel(:, 2), 'Color', '#EDB120');
xlabel('Time (s)');
ylabel('Velocity east (m/s)');
lgd = legend([p1, p2, p3], {'GNSS', 'Dead Reckoning', 'DR/GNSS Integrated Navigation'});
lgd.FontSize = 8;
title('All');
hold off
saveas(gcf, 'output/velocity_east.eps', 'epsc');


%% Heading graph 
figure('Renderer', 'painters', 'Position', [10 10 540 300])
plot(time, heading);
xlabel('Time (s)');
ylabel('Angle to north (°)');
lgd = legend([p1, p2, p3], {'GNSS', 'Dead Reckoning', 'DR/GNSS Integrated Navigation'});
lgd.FontSize = 8;
title('Gyro-magnetic heading');
saveas(gcf, 'output/heading.eps', 'epsc');



