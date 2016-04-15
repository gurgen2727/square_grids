#number of spanning trees of grid graph
load octave/matrix.in
format long;
d = det(matrix(1:length(matrix)-1,1:length(matrix)-1));
printf("Count of spanning trees of grid graph = %i\n", d);
#det(matrix(2:length(matrix),2:length(matrix)))
#det(matrix(1:length(matrix)-1,2:length(matrix)))
#det(matrix(2:length(matrix),1:length(matrix)-1))
nm = length(matrix);
m = 2;
n = 2;
for i = 3:nm 
	if( matrix(1,i) == -1 )
		n = i - 1;
		m = nm / n;
		break;
	endif
	if( matrix(i,1) == -1)
		m = i - 1;
		n = nm / m;
		break;
	endif
endfor

#print graphics of characteristic functions for hamiltonian loops
figure;
hold on;
grid on;
ylabel("Fc(k)");
tit = strcat("n=", num2str(n));
tit = strcat(tit,", m=");
tit = strcat(tit, num2str(m));
title(tit);
load octave/functions.in
[m n] = size(functions);
axis([1 n]);
color = [" - red,";" - green,";" - blue,";" - magenta,";" - cyan,";" - brown,"];
xlab = "k\n";
if( m < 8 )
	for i = 1:m
		xlab = strcat(xlab, strcat(" loop", num2str(i)));
		xlab = strcat(xlab, color(i,:));
	endfor
else 
	str = strcat("number of graphics=", num2str(m));
	xlab = strcat(xlab, str);
endif
xlabel(xlab);

col = "rgbmcb"; # red, green, blue, magenta, cyan, brown
for i = 1:m
	color = col( mod(i - 1, 6) + 1);
	opt = strcat(color, "o-.");
	plot(functions(i,:), opt);
endfor
pause;
