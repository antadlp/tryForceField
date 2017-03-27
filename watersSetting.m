function waterSetting(Org, volind, xs, ys, zs)

r = nthroot(3*volind/(4*pi), 3);

volTotal = xs*ys*zs;

a = 0;
b = 360;

xd = floor(xs/r);
yd = floor(ys/r);
zd = floor(zs/r);

numA = 3*xd*yd*zd;

file = fopen('waterSetting.xyz','w');
fprintf(file, '%d', numA);
fprintf(file, '\n\tmolec');
fprintf(file, '\n');

file2 = fopen('settingGeometry.dat','w');
file3 = fopen('settingConstantes.dat','w');

cteOH = 5.0;
cteHH = 10.0;

l=0;
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

         fprintf(file, 'O  %f  %f  %f\n', O(1), O(2), O(3));
         fprintf(file, 'H  %f  %f  %f\n', H1(1), H1(2), H1(3));
         fprintf(file, 'H  %f  %f  %f\n', H2(1), H2(2), H2(3));

        
         %Oxygen is bound only with two hidrogens
         
         fprintf(file2, '%d   %d',  ...
         l + 2, ...
         l + 3);
         fprintf(file2, '\n');
 

         fprintf(file3, '%d   %d',  ...
         cteOH, ...
         cteOH);
         fprintf(file3, '\n');
 

         %Hidrogens bound  with one oxygen and a rigid union
         %withe the other hidrogen
         % first %d : number of molecules to bound
         % second %d : oxigen to bound
         % third %d : hidrogen to bound 

         fprintf(file2, '%d   %d',  ...
         l + 1, ...
         l + 3);
         fprintf(file2, '\n');
 
         fprintf(file3, '%d   %d',  ...
         cteOH, ...
         cteHH);
         fprintf(file3, '\n');
 


         fprintf(file2, '%d   %d',  ...
         l + 1, ...
         l + 2);
         fprintf(file2, '\n');

         fprintf(file3, '%d   %d',  ...
         cteOH, ...
         cteHH);
         fprintf(file3, '\n');


         
         l = l + 3 ;




      end
   end
end

fclose(file);
fclose(file2);



         
               





