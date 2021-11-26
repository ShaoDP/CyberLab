close all;
figure
plot(out.y(:,1), out.x(:,1),out.y(:,2), out.x(:,2));axis equal;
legend("model 1", "model 2");