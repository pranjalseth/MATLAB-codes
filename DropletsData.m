clc
clear all
close all

%% reading the images
A=imread('DSC_0233.jpg');
figure(1),imshow(A);
B = (rgb2gray(A));
B=imsharpen (B,'Radius',1,'Amount',2);
B=adapthisteq(B);
figure(2),imshow(B);
C = im2bw(B,0.5);
figure(3),imshow(C);
% check here
B = bwareaopen(C,100);
figure(4),imshow(B);

%%
n1 = 233; %start image sequence
n2 = 301; %end image sequence
y2 = 3350; %estimate of the base, CHECK!

CR = zeros(n2-n1+1,1); %matrices to store values
H = zeros(n2-n1+1,1);
Ltheta = zeros(n2-n1+1,1);
Rtheta = zeros(n2-n1+1,1);
Dmax = zeros(n2-n1+1,1);

%% loop to read the entire run
for ni = n1:n2
	fprintf('%d',ni);
	if ni<10
		filename=sprintf('DSC_000%d.jpg',ni);
		A=imread(filename);
	else if ni<100
		filename=sprintf('DSC_00%d.jpg',ni);
		A=imread(filename);
	else if ni<1000
		filename=sprintf('DSC_0%d.jpg',ni);
		A=imread(filename);
	else if ni<10000
		filename=sprintf('DSC_%d.jpg',ni);
		A=imread(filename);
			% else
			% filename=sprintf('ImgA0%d.jpg',ni); 
			% A=imread(filename);
	end
	end
	end
	end
			%B = edge(rgb2gray(A),'sobel'); 

	%%
	B = im2bw(A,0.4); %MASTERCONTROL
	%figure(5), imshow(B);
	[m,n] = size(B);
	for i = 1:n
		for j = y2:m
			B(j,i) = 1;
		end
	end

	%figure(6),imshow(B);
	B=bwareaopen(B,10000);
	B = imcomplement(B);
	%figure(7),imshow(B);
	B = bwareaopen(B,10000);
	%figure(8),imshow(B);
	B = imfill(B,'holes');
	figure(9),imshow(B); %only leaves the droplet in white
	%white=1
	%black=0

	%% contact dia
	flag = 0;
	for i = 4000:-1:1
		for j = 1:6000
			if B(i,j) == 1;
				r1 = j;
				yr1 = i;
				flag = 1;
				break;
			end
		end
		if flag == 1;
			break;
		end
	end
	flag = 0;
	for i = 4000:-1:1
		for j = 6000:-1:1
			if B(i,j) == 1;
				r2 = j;
				yr2 = i;
				flag = 1;
				break;
			end
		end
		if flag == 1;
			break;
		end
	end
	r = r2-r1;
	CR(ni) = r;

	%% Height
	flag = 0;
	for i = 1:4000
		for j = 1:6000
			if B(i,j) == 1;
				h1 = i;
				xh1 = j;
				flag = 1;
				break;
			end
		end
		if flag == 1;
			break;
		end
	end
	flag = 0;
	for i = 4000:-1:1
		for j = 1:6000
			if B(i,j) == 1;
				h2 = i;
				xh2 = j;
				flag = 1;
				break;
			end
		end
		if flag == 1;
			break;
		end
	end
	h = h2-h1;
	H(ni) = h;

	%% Left CA
	anglefind1=0;
	start1 = r1;
	moveL = start1-10;
	moveR = start1+10;
	for i = yr1:-1:h1
		if B(i, moveL)==1
			l1=i;
			anglefind1=1;
			break;
		end
	end
	if anglefind1==0
		for i = yr1:-1:h1
			if B(i, moveR)==0
				l1=i;
				anglefind1=2;
				break;
			end
		end
	end
	if anglefind1==0 %90de
		ltheta = 90;
	end
	if anglefind1==1 %obtuse
		ltheta1 = atan((l1-h2)/10);
		ltheta = ltheta1*180/3.14;
		ltheta = 180+ltheta;
	end
	if anglefind1==2 %acute
		ltheta1 = atan((h2-l1)/10);
		ltheta = ltheta1*180/3.14;
	end
	Ltheta(ni) = ltheta;

	%% Right CA
	anglefind2=0;
	start2 = r2;
	moveL = start2-10;
	moveR = start2+10;
	for i = yr2:-1:h1
		if B(i, moveR)==1
			l2=i;
			anglefind2=1;
			break;
		end
	end
	if anglefind2==0
		for i = yr2:-1:h1
			if B(i, moveL)==0
				l2=i;
				anglefind2=2;
				break;
			end
		end
	end
	if anglefind2==0 %90de
		rtheta = 90;
	end
	if anglefind2==1 %obtuse
		rtheta1 = atan((l2-h2)/10);
		rtheta = rtheta1*180/3.14;
		rtheta = 180+rtheta;
	end
	if anglefind2==2 %acute
		rtheta1 = atan((h2-l2)/10);
		rtheta = rtheta1*180/3.14;
	end
	Rtheta(ni) = rtheta;

	%% maxdia
	if rtheta>90 && ltheta>90
		flag = 0;
		for j = 1:6000
			for i = 1:4000
				if B(i,j) == 1;
					dmax1 = j;
					ydmax1 = i;
					flag = 1;
					break;
				end
			end
			if flag == 1;
				break;
			end
		end
		flag = 0;
		for j = 6000:-1:1
			for i = 1:4000
				if B(i,j) == 1;
					dmax2 = j;
					ydamx2 = i;
					flag = 1;
					break;
				end
			end
			if flag == 1;
				break;
			end
		end
		dmax = dmax2-dmax1;
		Dmax(ni) = dmax;
	end
end
