function out=surff(n_dim,fskip,n1,n2)
%n_dim=2;   %1=density, 2=concentration
fnm='out'; %file sequence names
%fskip=100;
%% 
sizex=512; %array size in 2nd index; x-left/right in code
sizey=512; %array size in 1st index; y-up/down in code
for i=1:sizex
    x(i)=i-1; %physical size in x(left/right in code) 
end; 
for j=1:sizey
    y(j)=j-1; %physical size in y (up/down in code)
end;
%%
ymin=0;
ymax=sizey-1; %relevant physical size in up/down dir
xmin=0;
xmax=sizex-1; %rellavent physical size in left/right dir
if(n_dim==2)
 zmin=-0.165;
 zmax=-0.075;
else
  zmin=-1;
  zmax=1;
end;

%%
for f_num=n1:fskip:n2

figure %open figure window`

%"fch" string is the file being opened
fch=[fnm,'_',num2str(f_num)]
fich=load(fch);
%assign file to array "den"
den=fich(:,n_dim);

%plot data
rho=reshape(den, [sizey sizex]);
surf(x,y,rho,rho);
view(0,90);
caxis([zmin zmax]);
hold on;  
colorbar   
lighting phong;
shading interp;
%material shiny
%no grid;
%light('Position',[0 0 1]);

%determine which field is being plotted and thus labeled
if(n_dim==1)
 title(['Density. Time=',num2str(f_num)]);
end
if(n_dim==2)
 title(['Concentration. Time=',num2str(f_num)]); 
end

%set plot dimentions and window aspect ratio
axis([xmin xmax ymin ymax zmin zmax])
set(gca,'PlotBoxAspectRatio',[sizex,sizey,sizey])
grid on
xlabel('x-axis');
ylabel('y-axis');

   %%%printing figure to jpg file%%%%%%
   if(n_dim==1)
   fchnp1=['density',num2str(f_num),'.jpg'];
   end
   if(n_dim==2)
   fchnp1=['concentration',num2str(f_num),'.jpg'];
   end
   print('-r72','-djpeg100',fchnp1);
   close;

end;
%%
