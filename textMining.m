clear all
clc
close all


fid = fopen('C:\Users\babsf\Desktop\558PA1data.txt');
data = textscan(fid,'%s');
fclose(fid);

Gene=[];
A = data{1};
%A{1};
%for loop to split all 181
for i = 1 : 181
ggg = strsplit(A{i},',');
%remove last element from all, last element is a string
ggg(501) = [];

%read them into an array 181 x 500
gen = cellfun(@str2num,ggg);
Gene = [Gene; gen];
end

%Select five attributes assigned.
xdatatemp = Gene(:,[1 100 200 300 450]);

xx = xdatatemp';

minMatr =[];
for i = 1 : 5
    temp= min(xx(i,:));
    minMatr = [minMatr temp];
end
m = min(minMatr);

maxMatr =[];
for i = 1 : 5
    temp= max(xx(i,:));
    maxMatr = [maxMatr temp];
end
M = max(maxMatr);
rangeMatr = (M - m)/4;
g1= [];
g2= [];
g3= [];
g4= [];
for c = 1:5
    for r = 1:181
        
        if xx(c, r) < rangeMatr
            g1 = [g1 xx(c, r)];
        elseif xx(c, r) < (rangeMatr *2)
            g2 = [g2 xx(c, r)];
        elseif xx(c, r) < (rangeMatr *3)
            g3 = [g3 xx(c, r)];
        else
            g4 = [g4 xx(c, r)];
        end   
    end
end


%equal frequency approach
gg1= [];
gg2= [];
gg3= [];
gg4= [];
for c = 1:5
    [sortxx,Ixx] = sort(xx(c,:));
    for r = 1:181
        %sort
        if r < 46
            gg1 = [gg1 sortxx(r)];
        elseif r < 91
            gg2 = [gg2 sortxx(r)];
        elseif r < 136
            gg3 = [gg3 sortxx(r)];
        else
            gg4 = [gg4 sortxx(r)];
        end   
    end
end

gggg=[gg1; gg2; gg3];
%============================

 Outvalue = [];
 Outindex = [];
for i = 1 : 5
    y= (xx(i,:));
    [sorted_y,I] = sort(y);
    tempvalue = mat2cell(sorted_y, 1, [45 45 45 46]);
    tempindex = mat2cell(I, 1, [45 45 45 46]);
%     tempvalue = cellfun(@str2num,tempvalue);
%     tempindex = cellfun(@str2num,tempindex);
    Outvalue = [Outvalue; tempvalue];
    Outindex = [Outindex; tempindex];
end

%=======================================================
%choose 10 genes
T = Gene';
T = T(:,[5 23 44 55 57 97 111 124 130 141]);
T = T';

%find euclidean distance
D = pdist(T);
EuclideanDist = squareform(D);
me = 10000;
Me = 0;
for c = 1:10
    for r = 1:10
        if EuclideanDist(c,r) > M 
            Me =  EuclideanDist(c,r);
        end 
        if EuclideanDist(c,r) < me && EuclideanDist(c,r) ~= 0
            me =  EuclideanDist(c,r);
        end 
    end
end

[IMe,JMe] = find(EuclideanDist==Me);
[Ie,Je] = find(EuclideanDist==me);
%find cosine similarity distance
D2 = pdist(T,'cosine');
cosineSim = squareform(D2);
mc = 10000;
Mc = 0;
for c = 1:10
    for r = 1:10
        if cosineSim(c,r) > Mc 
            Mc =  cosineSim(c,r);
        end 
        if cosineSim(c,r) < mc && cosineSim(c,r) ~= 0
            mc =  cosineSim(c,r);
        end 
    end
end

[IMc,JMc] = find(cosineSim==Mc);
[Ic,Jc] = find(cosineSim==mc);
%correlation
%
%choose another 10 genes
T2 = Gene';
T2 = T2(:,[3 11 40 53 60 70 73 120 138 151]);
T2 = T2';
R = corr(T2');


%Apply min-max normalization to all attributes in the data set and repeat parts 
%min-max normalization 

tt =[];
for x = 1 : 10
    temp = -1 + 2.*(T(x,:) - min(T(x,:)))./(max(T(x,:)) - min(T(x,:)));
    %temp = sum(minus(T(x,:),lbpFeatures).^2);
    tt = [tt; temp];
end
%find euclidean distance
d = pdist(tt);
EuclideanDist2 = squareform(d);

%find cosine similarity distance
dd2 = pdist(tt,'cosine');
cosineSim2 = squareform(dd2);


%correlation
%min-max normalization 
tt2 =[];
for x = 1 : 10
    temp = -1 + 2.*(T2(x,:) - min(T2(x,:)))./(max(T2(x,:)) - min(T2(x,:)));
    %temp = sum(minus(T(x,:),lbpFeatures).^2);
    tt2 = [tt2; temp];
end
R2 = corr(tt2');
mc = 10000;
Mc = 0;
for c = 1:10
    for r = 1:10
        if R(c,r) > Mc && R(c,r) ~= 1
            Mc =  R(c,r);
        end 
        if R(c,r) < mc 
            mc =  R(c,r);
        end 
    end
end

[IMcor,JMcor] = find(R==M);
[Imcor,Jmcor] = find(R==m);

figure('Name','Equal-width plot');
hold on
plot(g1, 'k.')
plot(g2, 'b.')
plot(g3, 'g.')
plot(g4, 'r.')
hold off

figure('Name','Equal frequency plot');
hold on
plot(gg1, 'k.')
plot(gg2, 'b.')
plot(gg3, 'g.')
plot(gg4, 'r.')
hold off

 figure('Name','Correlations');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(R(:,plotId)) ;
 end
 
 figure('Name','Correlation after Normalization');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(R2(:,plotId)) ;
 end
 
   figure('Name','Euclidean Distances');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(EuclideanDist(:,plotId),'DisplayName','EuclideanDist') ;
 end
 

 figure('Name','Euclidean Distances after Normalization');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(EuclideanDist2(:,plotId)) ;
 end
 
  figure('Name','cosine Similarity ');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(cosineSim(:,plotId),'DisplayName','EuclideanDist') ;
 end
 

 figure('Name','cosine Similarity after Normalization');
 for plotId = 1 : 10
    subplot(3, 4, plotId) ;
     bar(cosineSim2(:,plotId)) ;
 end