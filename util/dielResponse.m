% $Id:$
function f=dielResponse(filename1,filename2,efield,ignore_dir)
%filename1 -> no bias
%filename2 -> with bias

%some parameters
unit_cell=5.13 %to average over
fsize=20

% open files
file1 =fopen(filename1,'r');
if file1<0
  disp('cannot open file '),
  disp(filename1)
  return
end

file2 =fopen(filename2,'r');
if file2<0
  disp('cannot open file '),
  disp(filename2)
  return
end

% read origin and end
orig1=fscanf(file1,'%lf',3);
end1 =fscanf(file1,'%lf',3);
fgetl(file1);

orig2=fscanf(file2,'%lf',3);
end2 =fscanf(file2,'%lf',3);
fgetl(file2);

% lattice constants
L1(1)=end1(1)-orig1(1)
L1(2)=end1(2)-orig1(2)
L1(3)=end1(3)-orig1(3)

L2(1)=end2(1)-orig2(1)
L2(2)=end2(2)-orig2(2)
L2(3)=end2(3)-orig2(3)

% Read grid size
np1=fscanf(file1,'%d',3)
fgetl(file1);
np2=fscanf(file2,'%d',3)
fgetl(file2);

% compute h
h(1)=L1(1)/np1(1);
h(2)=L1(2)/np1(2);
h(3)=L1(3)/np1(3);


nn3=np1(3)*np1(2)*np1(1);

% Read function
data1=fscanf(file1,'%lf',nn3);
data2=fscanf(file2,'%lf',nn3);

%conversion Vtotal into Hartree
data1=data1*0.5;
data2=data2*0.5;

%compute difference between data sets
diff=zeros(1,nn3);
for i=1:nn3
  diff(i)=data2(i)-data1(i);
end


%compute average over ignore_dir
dir1=1;
dir2=2;
dir3=ignore_dir

inc1=np1(2)*np1(3);
inc2=np1(3);
inc3=1

if (ignore_dir<3)
  dir2=3;
  inc2=1;
  inc3=np1(3);
end
if (ignore_dir<2)
  dir1=2;
  inc1=np1(3);
  inc3=np1(2)*np1(3);
end

nn2=np1(dir1)*np1(dir2);
avgdiff=zeros(1,nn2);
avg1   =zeros(1,nn2);
avg2   =zeros(1,nn2);
for i2=0:np1(dir2)-1
  for i1=0:np1(dir1)-1
    i=1+i1*inc1+inc2*i2;
    j=1+i1*np1(dir2)+i2;
    avgdiff(j)=0.;
    avg1(j)=0.;
    avg2(j)=0.;
    for i3=0:np1(dir3)-1
       avgdiff(j)=avgdiff(j)+diff(i3*inc3+i);
       avg1(j)=avg1(j)+data1(i3*inc3+i);
       avg2(j)=avg2(j)+data2(i3*inc3+i);
    end
    avgdiff(j)=avgdiff(j)/np1(dir3);
    avg1(j)=avg1(j)/np1(dir3);
    avg2(j)=avg2(j)/np1(dir3);
  end
end


% compute average over unit cells
cellradius=round(0.5*unit_cell/h(dir1))
invn2cell=1./(double((2*cellradius+1)*(2*cellradius+1)));

avgdiffcell=zeros(1,nn2);
for i1=cellradius:np1(dir1)-1-cellradius
  for i2=cellradius:np1(dir2)-1-cellradius
    l=i1*np1(dir2)+i2+1;
    
    % loop over mesh in unit cell
    avgdiffcell(l) = 0.;
    for i=-cellradius:cellradius
      for j=-cellradius:cellradius
        %k=l+i*inc1+j*inc2;
        k=l+i*np1(dir2)+j;
        avgdiffcell(l) = avgdiffcell(l)+avgdiff(k);
      end
    end
    avgdiffcell(l) = avgdiffcell(l)*invn2cell;
  end
end

nfdpts=2; % 2 for 4th order
%compute d/dx (interior points only)
ddxcell=zeros(1,nn2);
offset=cellradius+1; %+1 to avoid cells including boundary
invh1=1./h(dir1);
for i1=offset+nfdpts:np1(dir1)-nfdpts-offset
  for i2=offset:np1(dir2)-1-offset
    i=i1*np1(dir2)+i2+1;
    
    %staggered 2nd order FD
    %ddx(i)=(avgdiff(i)-avgdiff(i-np1(2)));
    
    %4th order FD
    %ddx(i)=(8./12.)*(avgdiff(i+np1(2))-avgdiff(i-np1(2)));
    %ddx(i)=ddx(i)+(1./12.)*(avgdiff(i-2*np1(2))-avgdiff(i+2*np1(2)));
    
    %staggered 4th order FD
    d1=(9./8.)*(avgdiffcell(i)-avgdiffcell(i-np1(dir2)));
    d2=(1./24.)*(avgdiffcell(i-2*np1(dir2))-avgdiffcell(i+np1(dir2)));
    
    ddxcell(i)=(d1+d2)*invh1;
    
    %add efield if using Hartree
    %ddxcell(i)=ddxcell(i)+efield;
    
    ddxcell(i)=efield/ddxcell(i);    
  end
end

%compute d/dx (interior points only)
ddx=zeros(1,nn2);
offset=1; %+1 to avoid boundary
for ix=offset+nfdpts:np1(dir1)-nfdpts-offset
  for iy=offset:np1(dir2)-1-offset
    i=ix*np1(dir2)+iy+1;
    
    %staggered 2nd order FD
    %ddx(i)=(avgdiff(i)-avgdiff(i-np1(2)));
    
    %4th order FD
    %ddx(i)=(8./12.)*(avgdiff(i+np1(2))-avgdiff(i-np1(2)));
    %ddx(i)=ddx(i)+(1./12.)*(avgdiff(i-2*np1(2))-avgdiff(i+2*np1(2)));
    
    %staggered 4th order FD
    ddx(i)=(9./8.)*(avgdiff(i)-avgdiff(i-np1(dir2)));
    ddx(i)=ddx(i)+(1./24.)*(avgdiff(i-2*np1(dir2))-avgdiff(i+np1(dir2)));
    
    ddx(i)=ddx(i)*invh1;

    %add efield if using Hartree
    %ddx(i)=ddx(i)+efield;

    ddx(i)=efield/ddx(i);
  end
end


% CONTOUR PLOTS in x,y plane
Xinit=orig1(dir1)
Xend =Xinit+L1(dir1)-h(dir1)

Yinit=orig1(dir2)
Yend =Yinit+L1(dir2)-h(dir2)

xcontour= linspace(Xinit,Xend,np1(dir1));
ycontour= linspace(Yinit,Yend,np1(dir2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxD=max(avgdiffcell)
minD=min(avgdiffcell)
v=linspace(minD,maxD,60);
%v=linspace(0.7,30,20)
W=reshape(avgdiffcell,[np1(dir2) np1(dir1)]); %matrix np1(1) x np1(2)

figure(8)
contourf(xcontour,ycontour,W,v)
title('Avg. over unit cells Diff. Vtotal)')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
axis equal
print('-depsc','-r100','picture8_dpi100')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maxD=min(10.,max(ddx))
%minD=max(min(ddx),0.5)
maxD=max(ddx)
minD=min(ddx)
%maxD=5.8
%minD=0.
v=linspace(minD,maxD,50);
v=linspace(0.,7.2,20);
W=reshape(ddx,[np1(dir2) np1(dir1)]);
%W=W';

figure(10)
contourf(xcontour,ycontour,W,v)
title('Efield/(d/dx average diff)')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxD=max(ddx)
minD=min(ddx)
v=linspace(0.7,1.3,20);
W=reshape(ddx,[np1(dir2) np1(dir1)]);

figure(11)
contourf(xcontour,ycontour,W,v)
title('blow up Efield/(d/dx average diff)')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxD=max(avgdiff)
minD=min(avgdiff)
%v=linspace(-L2(1)*efield*0.5,L2(1)*efield*0.5,50);
v=linspace(minD,maxD,50);
W=reshape(avgdiff,[np1(dir2) np1(dir1)]);

figure(12)
contourf(xcontour,ycontour,W,v)
title('z average diff')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
print('-depsc','-r100','picture_dpi100')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxD=max(avg1)
minD=min(avg1)
v=linspace(minD,maxD,20);
%v=linspace(-0.001,0.001,20);
W=reshape(avg1,[np1(dir2) np1(dir1)]);

figure(13)
contourf(xcontour,ycontour,W,v)
title('z average data1')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
print('-depsc','-r100','pictureData1_dpi100')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxD=max(avg2)
minD=min(avg2)
v=linspace(minD,maxD,20);
range=L2(dir1)*efield
%v=linspace(-L2(1)*efield*0.5,L2(1)*efield*0.5,20)
W=reshape(avg2,[np1(dir2) np1(dir1)]);

figure(14)
contourf(xcontour,ycontour,W,v)
title('z average data2')
axis([Xinit Xend Yinit Yend])
colorbar
xlabel('x')
ylabel('y')
print('-depsc','-r100','pictureData2_dpi100')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=reshape(ddxcell,[np1(dir2) np1(dir1)]);
nbuffer1=cellradius+1        % rows of matrix (y)
nbuffer2=cellradius+nfdpts+1 % columns (x)
size(W)
Z=W(1+nbuffer1:np1(dir2)-nbuffer1,1+nbuffer2:np1(dir1)-nbuffer2);

maxD=max(Z(:))
minD=min(Z(:))

xi=Xinit+nbuffer1*h(dir1)
xe=Xend -nbuffer1*h(dir1)
yi=Yinit+nbuffer2*h(dir2)
ye=Yend -nbuffer2*h(dir2)
[r,c]=size(Z)
xcontour= linspace(xi,xe,c);
ycontour= linspace(yi,ye,r);

figure(9)
ncontours = 20;

vmin=0.;
vmax=6.5;
v=linspace(vmin,vmax,ncontours)

contourf(xcontour,ycontour,Z,v)
set(gca,'fontsize',fsize)

%title('Epsilon')

caxis([vmin vmax])

colormap(jet)
axis([xi xe yi ye])
colorbar('FontSize',fsize)

xlabel('x','fontsize',fsize)
ylabel('y','fontsize',fsize)
axis equal

print('-depsc','-r100','picture_dpi100') 
clear

