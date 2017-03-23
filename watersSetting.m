function waterSetting(Org, volind, xs, ys, zs)

r = nthroot(3*volind/(4*pi), 3);

volTotal = xspace*yspace*zspace;

a = 0;
b = pi;

xd = floor(xs/r);
yd = floor(ys/r);
zd = floor(zs/r);

for i=1:xd
   for j=1:yd
      for k=1:zd
         
         rx = a + (b-a).*rand(1, 1);
         ry = a + (b-a).*rand(1, 1);
         rz = a + (b-a).*rand(1, 1);

         xo = Org(1) + (i-1)*r;
         yo = Org(2) + (j-1)*r;
         zo = Org(3) + (k-1)*r;

         Or = [xo, yo, zo];

         [O H1 H2] = singleWaterTemplate(Or, rx, ry, rz);







         
               





