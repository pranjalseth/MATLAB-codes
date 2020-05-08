clc
clear all
close all

%4-4 20,rbc
%% reading the images
A=imread('4-4-600-rbcs20um-c00001.tif');
figure(1),imshow(A);

B = im2bw(A,0.38);   %master-control
figure(2),imshow(B);

C = imcomplement(B);
figure(3),imshow(C);
C = bwareaopen(C,60);
figure(4),imshow(C);
% C = imfill(C,'holes');
% figure(5),imshow(C);

%white=1
%black=0

%%
BottomL = [75,650];
BottomR = [860,650];
TopL = [75,150];
TopR = [860,150];
Xmin = 75;
Xmax = 860;
Ymax = 650;
Ymin = 150;

%%

L = Xmax-Xmin;
H = Ymax-Ymin;


D = imcrop(C,[Xmin Ymin L H]);
imshow(D);

E = bwareaopen(D,50);
figure(6),imshow(E);

stats = regionprops(E,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
p = size(stats);
q = p(1);

for i = 1:1:q
    center = stats(i).Centroid;
    Smax = stats(i).MajorAxisLength;
    Smin = stats(i).MinorAxisLength;
    X = 0.93*center(1);
    Y = center(2);
    Y = 0.93*(500-Y);
    
    if Smax<=15 && Smin>=5
        figure(9), scatter(X,Y,'b');
        axis([0 720 0 480]);
        hold on;
    end
    
    if Smax<=40 && Smin>=10
        figure(9), scatter(X,Y,'r')
        axis([0 720 0 480])
        hold on
    end
    
end


% stats = regionprops(E,'Centroid');
% center = stats.Centroid;


%   ALL THIS WAS A TEST!
%   NOW SHOWTIME!


%%
%%
%%
%%


%SHOWTIME!!!

n1 = 1; %start image sequence
n2 = 150; %24848; %end image sequence

centers = zeros(n2-n1+1,2); %matrix to store values

%%
% to import image sequence

for ni = n1:1:n2
     fprintf('%d',ni);
        if ni<10
            filename=sprintf('4-4-600-rbcs20um-c0000%d.tif',ni);
            A=imread(filename);
        else if ni<100
            filename=sprintf('4-4-600-rbcs20um-c000%d.tif',ni);
            A=imread(filename);
        else if ni<1000
            filename=sprintf('4-4-600-rbcs20um-c00%d.tif',ni);
            A=imread(filename);
        else if ni<10000
            filename=sprintf('4-4-600-rbcs20um-c0%d.tif',ni);
            A=imread(filename);
        else if ni<100000
            filename=sprintf('4-4-600-rbcs20um-c%d.tif',ni);
            A=imread(filename);
            
        end
        end
        end 
        end
        end
        
        

B = im2bw(A,0.35);
C = imcomplement(B);
C = bwareaopen(C,60);
% C = imfill(C,'holes');
%white=1
%black=0

D = imcrop(C,[Xmin Ymin L H]);
E = bwareaopen(D,50);
figure(8),imshow(E);

stats = regionprops(E,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
p = size(stats);
q = p(1);

for i = 1:1:q
    center = stats(i).Centroid;
    Smax = stats(i).MajorAxisLength;
    Smin = stats(i).MinorAxisLength;
    X = 0.93*center(1);
    Y = center(2);
    Y = 0.93*(500-Y);
    
    if Smax<=15 && Smin>=5
        figure(9), scatter(X,Y,'b')
        axis([0 720 0 480])
        hold on
    end
    
    if Smax<=40 && Smin>=10
        figure(9), scatter(X,Y,'r')
        axis([0 720 0 480])
        hold on
    end
    
end

% stats = regionprops(E,'Centroid');
% center = stats.Centroid;
% 
% centers(ni,1) = center(1);
% centers(ni,2) = center(2);

end

% X = 0.93*centers(:,1);
% Y = centers(:,2);
% Y = 0.93*(515-Y);
% scatter(X,Y);
% axis([0 720 0 480])
