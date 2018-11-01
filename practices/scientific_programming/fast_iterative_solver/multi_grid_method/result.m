
err = importdata('p2_6/error1.txt');
err2 = importdata('p2_6/error.txt');

r_time = importdata('p2_6/r_time1.txt');
r_time2 = importdata('p2_6/r_time.txt');

x = linspace(0,length(err)-1,length(err));

semilogy(x,err,'*-','markersize',10);
hold on;
semilogy(x,err2,'d-','markersize',10);
legend('\nu1=\nu2=1','\nu1=2, \nu2=1');
title('Error graph (n=4)');
xlabel('The number of iteration');
ylabel('Relative error');
grid;
hold off;

