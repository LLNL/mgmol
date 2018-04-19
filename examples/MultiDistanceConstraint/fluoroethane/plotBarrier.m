function f=plotBarrier(filename)

% open data file
file1 =fopen(filename,'r');
if file1<0
  disp('cannot open file '),
  disp(filename)
  return
end

%skip first line (comments)
tline = fgetl(file1);

%read three floats per line and skip two strings
C = textscan(file1,'%f %f %f %*s %*s')

M = cell2mat(C);
M = sortrows(M,1)

%C=mat2cell(M);
C={M};
c1=M(:,1);
c2=M(:,2);
c3=M(:,3);

ha2kcalpmol=627.51057;
bohr2ang=0.529177;

%
ref_energy=c2(1);

c2=c2-ref_energy;

%convert energies and forces to kcal/mol
c2=c2*ha2kcalpmol;
c3=c3*ha2kcalpmol;

%convert distances to Ang
c1=c1*bohr2ang;
c3=c3/bohr2ang;

sp1=size(c1);
m1=sp1(1);

n=2*3;

f1=zeros(n*(m1-1)+1,2);
f1(1,1)=c1(1);
f1(1,2)=0.;
val=0.;
for i=1:m1-1
  dx=c1(i+1)-c1(i);
  if dx>0.
    c=val;
    b=c3(i);
    a=(c3(i+1)-b)/(2.*dx);
    for j=1:n
      d=j*dx/n;
      x=c1(i)+d;
      val=a*d*d+b*d+c;
      f1(n*(i-1)+j+1,1)=x;
      f1(n*(i-1)+j+1,2)=val;
    end
  end
end


f1

plot(f1(:,1),f1(:,2))
%text(f1(m1*n/2,1),f1(m1*n/2,2),'transition state')
text(f1(1,1),f1(1,2),'product')
text(f1(n*(m1-1),1),f1(n*(m1-1),2),'reactants')

hold on

plot(c1(:),c2(:),'r+')

hold on

xlabel('reaction coord. [Ang]')
ylabel('energy [kcal/mol]')

print('-depsc','-r100','plotBarrier_dpi100') 
clear

