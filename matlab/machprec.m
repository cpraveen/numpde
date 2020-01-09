h=1.0/2.0;
while 1+h > 1
   h = h/2;
end
h = 2*h;
fprintf(1,'Machine precision = %e\n', h)
