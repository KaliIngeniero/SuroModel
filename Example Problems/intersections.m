function Dist=intersections(vars)
% Dist=intersections(vars)
%
% Calculates level of intersections in a novel vibration isolating
% structure
%
% Inputs: 
%   vars - n x 18,36,54,... vector of joint positions with elements in range [0,1] 
%
% Outputs: 
%   Dist - scalar total intereference of beam elements
%
% Copyright 2007 A I J Forrester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

warning off MATLAB:divideByZero
% Transpose vars if necessary
if size(vars,1)~=1
	vars=vars';
end
Length=0.3;
Range=0.125;
EndTol=0.01;
Goal=0.0075;
%%% Build matrix of normalised boom node positions %%%
for i=1:length(vars)/3
NormVars(i,1:3)=vars((i-1)*3+1:(i-1)*3+3);
end
%%% Build matrix of non normalised boom node positions %%%
for k=1:size(NormVars,1)/3
	for j=1:3
		if j==1
			if round(k/2)==k/2
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+EndTol;
			else
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+Length-Range-EndTol;
			end
			Vars((k-1)*3+j,2)=NormVars((k-1)*3+j,2)*Range-Range/2;
			Vars((k-1)*3+j,3)=NormVars((k-1)*3+j,3)*Range+Length-Range/2;
		elseif j==2
			if round(k/2)==k/2
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+EndTol;
			else
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+Length-Range-EndTol;
			end
			Vars((k-1)*3+j,2)=NormVars((k-1)*3+j,2)*Range-Range/2;
			Vars((k-1)*3+j,3)=NormVars((k-1)*3+j,3)*Range-Range/2;
		else
			if round(k/2)==k/2
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+EndTol;
			else
				Vars((k-1)*3+j,1)=NormVars((k-1)*3+j,1)*Range+Length-Range-EndTol;
			end
			Vars((k-1)*3+j,2)=-1*(NormVars((k-1)*3+j,2)*-1*Range+(Length^2-(Length/2)^2)^0.5+Range/2);
			Vars((k-1)*3+j,3)=NormVars((k-1)*3+j,3)*Range+Length/2-Range/2;
		end
	end
end

%%% Fixed base coordinates %%%
Base=[0 0 Length
0 0 0
0 -(Length^2-(Length/2)^2)^0.5 Length/2];
%%% Fixed mount coordinates %%%
Mount=[Length 0.0 Length
Length 0.0 0.0
Length -(Length^2-(Length/2)^2)^0.5 Length/2];
Beam=[Base;Vars;Mount];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find ball end positions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StartBeam=Beam;
%%% for all nodes calculate ball joint angle %%
for i=1:size(Beam,1)
	if i<4 | i>size(Beam,1)-3
		if (i+2)/3==round((i+2)/3)
			AlphaVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:))+(Beam(i+2,:) - Beam(i,:))./norm(Beam(i+2,:) - Beam(i,:)))./2;	
		end
		if (i+1)/3==round((i+1)/3)
			AlphaVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:))+(Beam(i-1,:) - Beam(i,:))./norm(Beam(i-1,:) - Beam(i,:)))./2;
		end
		if (i)/3==round((i)/3)
			AlphaVector(i,:)=((Beam(i-1,:) - Beam(i,:))./norm(Beam(i-1,:) - Beam(i,:))+(Beam(i-2,:) - Beam(i,:))./norm(Beam(i-2,:) - Beam(i,:)))./2;
		end
		if i<4
			BETA(i)=pi/2;
			%JointVector(i,:)=[1 0 0];
		elseif i>size(Beam,1)-3
			%JointVector(i,:)=[-1 0 0];
			BETA(i)=-pi/2;
		end
	else
		if (i-1)/3==round((i-1)/3)
			JointVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:)) ...
			+ (Beam(i+2,:) - Beam(i,:))./norm(Beam(i+2,:) - Beam(i,:)) ...
			+ (Beam(i-3,:) - Beam(i,:))./norm(Beam(i-3,:) - Beam(i,:)) ...
			+ (Beam(i+3,:) - Beam(i,:))./norm( Beam(i+3,:) - Beam(i,:))...
			+ (Beam(i-1,:) - Beam(i,:))./ norm(Beam(i-1,:) - Beam(i,:))...
			+ (Beam(i+4,:) - Beam(i,:))./ norm(Beam(i+4,:) - Beam(i,:)))./6; 
			AlphaVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:))+(Beam(i+2,:) - Beam(i,:))./norm(Beam(i+2,:) - Beam(i,:)))./2;
		end
		if (i-2)/3==round((i-2)/3)
			JointVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:)) ...
			+( Beam(i-1,:)- Beam(i,:))./ norm(Beam(i-1,:)- Beam(i,:)) ...
			+ (Beam(i-3,:)- Beam(i,:))./ norm(Beam(i-3,:)- Beam(i,:)) ...
			+( Beam(i+3,:)- Beam(i,:))./ norm(Beam(i+3,:)- Beam(i,:)) ...
			+( Beam(i+4,:)- Beam(i,:))./ norm(Beam(i+4,:)- Beam(i,:)) ...
			+( Beam(i-4,:)- Beam(i,:))./ norm(Beam(i-4,:)- Beam(i,:)))./6;
			AlphaVector(i,:)=((Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:))+(Beam(i-1,:) - Beam(i,:))./norm(Beam(i-1,:) - Beam(i,:)))./2;
		end
		if i/3==round(i/3)
			JointVector(i,:)=((Beam(i-1,:)- Beam(i,:))./ norm(Beam(i-1,:)- Beam(i,:)) ... 
			+ (Beam(i-2,:) - Beam(i,:))./norm(Beam(i-2,:) - Beam(i,:)) ...
			+ (Beam(i-3,:) - Beam(i,:))./norm(Beam(i-3,:) - Beam(i,:)) ...
			+ (Beam(i+3,:) - Beam(i,:))./norm(Beam(i+3,:) - Beam(i,:)) ...
			+ (Beam(i-4,:) - Beam(i,:))./norm(Beam(i-4,:) - Beam(i,:)) ...
			+ (Beam(i+1,:) - Beam(i,:))./norm(Beam(i+1,:) - Beam(i,:)))./6;
			AlphaVector(i,:)=((Beam(i-1,:) - Beam(i,:))./norm(Beam(i-1,:) - Beam(i,:))+(Beam(i-2,:) - Beam(i,:))./norm(Beam(i-2,:) - Beam(i,:)))./2;
		end	
	end			
end

JointVector=[1 0 0; 1 0 0; 1 0 0; JointVector(4:end,:); -1 0 0; -1 0 0; -1 0 0];
[Gamma,BETA,DummyR]=cart2sph(JointVector(:,1),JointVector(:,2),JointVector(:,3));
[psi,delta,R]=cart2sph(AlphaVector(:,1),AlphaVector(:,2),AlphaVector(:,3));
BETA=-BETA;
Gamma=-Gamma;
JointVector;
AlphaVector;
%%% Calculate ball centres %%%
for i=1:size(Beam,1)
	THETA=psi(i);
	THETADASH=THETA+Gamma(i);
	A=sqrt(AlphaVector(i,1)^2+AlphaVector(i,2)^2)*cos(THETADASH);
	ALPHA=atan(AlphaVector(i,3)/A);
	DELTA=pi/2-ALPHA+BETA(i);
	DELTA=pi/2-ALPHA-BETA(i);
	D=sqrt(AlphaVector(i,3)^2+A^2);
	C=D*cos(DELTA);
	E=sqrt(AlphaVector(i,1)^2+AlphaVector(i,2)^2)*sin(THETADASH);
	ETA(i)=atan(E/C);
	deltadash(i)=ETA(i);
if i==1
deltadash(i)=-5*pi/6;
elseif i==2
deltadash(i)=-pi/6;
elseif i==3
deltadash(i)=pi/2;
elseif (i==4 | i==10) & Beam(i,1)==Beam(i+1,1) & Beam(i,1)==Beam(i+2,1)
deltadash(i)=pi+ETA(i);
elseif (i==5 | i==11) & Beam(i,1)==Beam(i+1,1) & Beam(i,1)==Beam(i-1,1)
deltadash(i)=ETA(i);
elseif (i==6 | i==12) & Beam(i,1)==Beam(i-1,1) & Beam(i,1)==Beam(i-2,1)
deltadash(i)=ETA(i);
elseif (i==7 | i==13) & Beam(i,1)==Beam(i+1,1) & Beam(i,1)==Beam(i+2,1)
deltadash(i)=pi+ETA(i);
elseif (i==8 | i==15) & Beam(i,1)==Beam(i+1,1) & Beam(i,1)==Beam(i-1,1)
deltadash(i)=ETA(i);
elseif (i==9 | i==16) & Beam(i,1)==Beam(i-1,1) & Beam(i,1)==Beam(i-2,1)
deltadash(i)=pi+ETA(i);
elseif i==size(Beam,1)-2
deltadash(i)=5*pi/6;
elseif i==size(Beam,1)-1
deltadash(i)=pi/6;
elseif i==size(Beam,1)
deltadash(i)=-pi/2;
end




if JointVector(i,1)>0 & JointVector(i,2)>0 & JointVector(i,3)>0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	if ETA(i)<0
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	if ETA(i)<0
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	end
elseif  JointVector(i,1)>0 & JointVector(i,2)<0 & JointVector(i,3)>0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0	
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	eltadash=ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	end
elseif  JointVector(i,1)<0 & JointVector(i,2)>0 & JointVector(i,3)>0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	if ETA(i)<0
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	if ETA(i)<0
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	end
elseif  JointVector(i,1)<0 & JointVector(i,2)<0 & JointVector(i,3)>0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	end
elseif  JointVector(i,1)>0 & JointVector(i,2)>0 & JointVector(i,3)<0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	if ETA(i)<0
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	if ETA(i)<0
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	end
	
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	end
elseif  JointVector(i,1)>0 & JointVector(i,2)<0 & JointVector(i,3)<0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=+ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	end
elseif  JointVector(i,1)<0 & JointVector(i,2)>0 & JointVector(i,3)<0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0 %changed aijf 05/04/06
	if ETA(i)<0
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	if ETA(i)<0
	deltadash(i)=ETA(i); %changed aijf 13/02/06
	else
	deltadash(i)=pi+ETA(i); %changed aijf 13/02/06
	end
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	end
elseif  JointVector(i,1)<0 & JointVector(i,2)<0 & JointVector(i,3)<0
	if AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)>0
	deltadash(i)=-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)>0
	deltadash(i)=ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)>0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	if ETA(i)<0
		deltadash(i)=pi+ETA(i); %changed aijf 09/10/07
	else
		deltadash(i)=pi-ETA(i); %changed aijf 09/10/07
	end
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)>0 & AlphaVector(i,3)<0
	deltadash(i)=pi-ETA(i);
	elseif AlphaVector(i,1)<0 & AlphaVector(i,2)<0 & AlphaVector(i,3)<0
	deltadash(i)=pi+ETA(i);
	end
end

	for j=1:6
		TiltedBall(j,1,i)=0.006*sin((pi/2)-pi/6-deltadash(i)-(j-1)*2*pi/6)*cos(pi/2-BETA(i));
		TiltedBall(j,2,i)=0.006*cos((pi/2)-pi/6-deltadash(i)-(j-1)*2*pi/6);
		TiltedBall(j,3,i)=0.006*sin((pi/2)-pi/6-deltadash(i)-(j-1)*2*pi/6)*sin(pi/2-BETA(i));
		aa(j,i)=TiltedBall(j,1,i)/cos(Gamma(i));
		bb(j,i)=aa(j,i)*sin(Gamma(i));
		RotatedBall(j,1,i)=(TiltedBall(j,2,i)-bb(j,i))*sin(Gamma(i))+aa(j,i);
		RotatedBall(j,2,i)=(TiltedBall(j,2,i)-bb(j,i))*cos(Gamma(i));
		RotatedBall(j,3,i)=TiltedBall(j,3,i);
	end
	BallGamma(:,1,i)=RotatedBall(:,1,i)+Beam(i,1);
	BallGamma(:,2,i)=RotatedBall(:,2,i)+Beam(i,2);
	BallGamma(:,3,i)=RotatedBall(:,3,i)+Beam(i,3);

end

BallEnd=BallGamma;	
%%% Make beam database %%%
CBeam=Beam;
CBeam(:,2)=CBeam(:,2)+abs((Beam(1,3)/2)*tan(30*pi/180));
CBeam(:,3)=CBeam(:,3)-Beam(1,3)/2;
Angle=cart2pol(CBeam(:,2),CBeam(:,3));
%vector of empty ball slots for "vertical" beams
VacantSlot=ones(size(Beam,1),1).*3;
B=0;
V=0;
%Vertical beam intersection elimination
for J=1:3
	B1=3;
	B2=4;
	I=1;

	for I=1:5
		for i=1:(size(Beam,1)-3)/3
			B=B+1;
			if I==1
				% 4-3-4-3 
				TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(3,:,(i-1)*3+J+3) BallEnd(4,:,(i-1)*3+J)];
			% 3-4-3-4 
			elseif I==2
				TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(4,:,(i-1)*3+J+3) BallEnd(3,:,(i-1)*3+J)];
			% 4-4-3-3
			elseif I==3
				if i==1 | i==3 | i==5 | i==7 | i==9| i==11 
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(4,:,(i-1)*3+J+3) BallEnd(4,:,(i-1)*3+J)];
				else
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(3,:,(i-1)*3+J+3) BallEnd(3,:,(i-1)*3+J)];
				end
			elseif I==4
			% 3-3-4-4
				if i==1 | i==3 | i==5 | i==7 | i==9| i==11 
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(3,:,(i-1)*3+J+3) BallEnd(3,:,(i-1)*3+J)];
				else
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(4,:,(i-1)*3+J+3) BallEnd(4,:,(i-1)*3+J)];
				end
			elseif I==5
			% 3-4-3-3-4-3
				if i==1
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(4,:,(i-1)*3+J+3) BallEnd(3,:,(i-1)*3+J)];				
				elseif i==2
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(3,:,(i-1)*3+J+3) BallEnd(3,:,(i-1)*3+J)];				
				elseif i==3
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(3,:,(i-1)*3+J+3) BallEnd(4,:,(i-1)*3+J)];			
				elseif i==4
					TentativeBeamDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[BallEnd(4,:,(i-1)*3+J+3) BallEnd(4,:,(i-1)*3+J)];	
				end
			end
			JointDataBase((J-1)*(size(Beam,1)-3)/3+i,:)=[Beam((i-1)*3+J+3,:) Beam((i-1)*3+J,:)];
		end
		for i=1:size(TentativeBeamDataBase,1)
			x1=TentativeBeamDataBase(i,1:3);
			x2=TentativeBeamDataBase(i,4:6);
			for j=1:size(TentativeBeamDataBase,1)
				x3=TentativeBeamDataBase(j,1:3);
				x4=TentativeBeamDataBase(j,4:6);
				a=x2-x1;
				b=x4-x3;
				c=x3-x1;
				if i~=j 
				% are lines planar ?
					if dot((x1-x3),cross((x2-x1),(x4-x3)))==0
						q1=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
						q2=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
						v=0;
					else
						%position of minimum distance
						u1=a/norm(a);
						u2=b/norm(b);
						u3=cross(u1,u2)/norm(cross(u1,u2));
						e=dot((x1-x3),u3)*u3;
						v=norm(dot(c,cross(a,b)))/norm(cross(a,b));
						x3dash=x3+e;
						x4dash=x4+e;
						a=x2-x1;
						b=x4dash-x3dash;
						c=x3dash-x1;
						q1=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
						e=dot((x3-x1),u3)*u3;
						x1dash=x1+e;
						x2dash=x2+e;				
						a=x2dash-x1dash;
						b=x4-x3;
						c=x3-x1dash;
						q2=x1dash+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
					end
					BalltoMinDist(1)=sqrt(sum((q1-x1).^2));
					BalltoMinDist(2)=sqrt(sum((q1-x2).^2));
					BalltoMinDist(3)=sqrt(sum((q1-x3).^2));
					BalltoMinDist(4)=sqrt(sum((q1-x4).^2));
					[a,b]=sort([x1(1) x2(1)]);
					[c,d]=sort([x3(1) x4(1)]);
					if q1(1)>a(1) & q1(1)<a(2) & q2(1)>c(1) & q2(1)<c(2)
						[a,b]=sort([x1(2) x2(2)]);
						[c,d]=sort([x3(2) x4(2)]);
						if q1(2)>a(1) & q1(2)<a(2) & q2(2)>c(1) & q2(2)<c(2)
							[a,b]=sort([x1(3) x2(3)]);
							[c,d]=sort([x3(3) x4(3)]);
							if q1(3)>a(1) & q1(3)<a(2) & q2(3)>c(1) & q2(3)<c(2)
						%%% Check if min distance is very close to joint and allow v>0.0055 if this is the case 
								if min(BalltoMinDist)<0.02
									if v < 0.0053
										V=V+0.0053-v;
									end
								else
									if v < Goal
										V=V+Goal-v;
									end
								end
							end
						end
					end			
				end
			end	
			if V==0; 
				I=6; 
			end
		end
		V=0;
	end
end
Vvertical=V;
BeamDataBase=TentativeBeamDataBase;

B=size(BeamDataBase,1);
for i=4:size(Beam,1)
	if round(i/3)==i/3
		B=B+1;
		if (i/3)==1 |(i/3)==3 |(i/3)==5 |(i/3)==7 |(i/3)==9|(i/3)==11 ...
		| ((i+1)/3)==1 |((i+1)/3)==3 |((i+1)/3)==5 |((i+1)/3)==7 |((i+1)/3)==9|((i+1)/3)==11 ...
		| ((i+2)/3)==1 |((i+2)/3)==3 |((i+2)/3)==5 |((i+2)/3)==7 |((i+2)/3)==9|((i+2)/3)==11 
			BeamDataBase(B,:)=[BallEnd(6,:,i) BallEnd(1,:,i-2)];
		else
			BeamDataBase(B,:)=[BallEnd(1,:,i) BallEnd(6,:,i-2)];
		end
		JointDataBase(B,:)=[Beam(i,:) Beam(i-2,:)];
	else
		B=B+1;
		if (i/3)==1 |(i/3)==3 |(i/3)==5 |(i/3)==7 |(i/3)==9|(i/3)==11 ...
		| ((i+1)/3)==1 |((i+1)/3)==3 |((i+1)/3)==5 |((i+1)/3)==7 |((i+1)/3)==9|((i+1)/3)==11 ...
		| ((i+2)/3)==1 |((i+2)/3)==3 |((i+2)/3)==5 |((i+2)/3)==7 |((i+2)/3)==9|((i+2)/3)==11 
		BeamDataBase(B,:)=[BallEnd(6,:,i) BallEnd(1,:,i+1)];
		else
		BeamDataBase(B,:)=[BallEnd(1,:,i) BallEnd(6,:,i+1)];
		end
		JointDataBase(B,:)=[Beam(i,:) Beam(i+1,:)];
	end
	if round((i-1)/3)==(i-1)/3
		B=B+1;
			if (i/3)==1 |(i/3)==3 |(i/3)==5 |(i/3)==7 |(i/3)==9|(i/3)==11 ...
	| ((i+1)/3)==1 |((i+1)/3)==3 |((i+1)/3)==5 |((i+1)/3)==7 |((i+1)/3)==9|((i+1)/3)==11 ...
	| ((i+2)/3)==1 |((i+2)/3)==3 |((i+2)/3)==5 |((i+2)/3)==7 |((i+2)/3)==9|((i+2)/3)==11 
		BeamDataBase(B,:)=[BallEnd(2,:,i) BallEnd(2,:,i-1)];
		else
		BeamDataBase(B,:)=[BallEnd(5,:,i) BallEnd(5,:,i-1)];
		end
		JointDataBase(B,:)=[Beam(i,:) Beam(i-1,:)];
	elseif round((i-2)/3)==(i-2)/3
		B=B+1;
		if (i/3)==1 |(i/3)==3 |(i/3)==5 |(i/3)==7 |(i/3)==9|(i/3)==11 ...
	| ((i+1)/3)==1 |((i+1)/3)==3 |((i+1)/3)==5 |((i+1)/3)==7 |((i+1)/3)==9|((i+1)/3)==11 ...
	| ((i+2)/3)==1 |((i+2)/3)==3 |((i+2)/3)==5 |((i+2)/3)==7 |((i+2)/3)==9|((i+2)/3)==11 
		BeamDataBase(B,:)=[BallEnd(2,:,i) BallEnd(2,:,i-4)];
		else
		BeamDataBase(B,:)=[BallEnd(5,:,i) BallEnd(5,:,i-4)];
		end
		JointDataBase(B,:)=[Beam(i,:) Beam(i-4,:)];
	else
		B=B+1;
		if (i/3)==1 |(i/3)==3 |(i/3)==5 |(i/3)==7 |(i/3)==9|(i/3)==11 ...
	| ((i+1)/3)==1 |((i+1)/3)==3 |((i+1)/3)==5 |((i+1)/3)==7 |((i+1)/3)==9|((i+1)/3)==11 ...
	| ((i+2)/3)==1 |((i+2)/3)==3 |((i+2)/3)==5 |((i+2)/3)==7 |((i+2)/3)==9|((i+2)/3)==11 
		BeamDataBase(B,:)=[BallEnd(2,:,i) BallEnd(2,:,i-4)];
		else
		BeamDataBase(B,:)=[BallEnd(5,:,i) BallEnd(5,:,i-4)];
		end
		JointDataBase(B,:)=[Beam(i,:) Beam(i-4,:)];
	
	end
end
%%% Concatonate with base and mount %%%	
	
% check if Joints coinside %
V=sum(sum(Beam))*0;
for i=1:size(Beam,1)
	for j=1:size(Beam,1)
		if i~=j
			JointSep=sqrt(sum((Beam(i,:)-Beam(j,:)).^2));
			if JointSep<0.0225
				V=V+abs(JointSep-0.0225);
			end
		end
	end
end
Vjoints=V-Vvertical;
%check if Joints lie on a beam
for i=1:size(BeamDataBase,1)
	for j=1:size(Beam,1)
		if JointDataBase(i,1)~=Beam(j,1) & JointDataBase(i,2)~=Beam(j,2) & JointDataBase(i,3)~=Beam(j,3)
			x1=BeamDataBase(i,1:3);
			x2=BeamDataBase(i,4:6);
			x0=Beam(j,:);
			d=norm(cross(x2-x1,x1-x0))/norm(x2-x1);
			if (x1<x2 & x0>x1 & x0<x2) | (x1>x2 & x0<x1 & x0>x2);
				if d < 0.01125+Goal/2
					V=V+0.01125+Goal-d;
				end
			end
		end
	end
end	
Vjoints2=V-Vvertical-Vjoints;

for i=1:size(BeamDataBase,1)
	x1=BeamDataBase(i,1:3);
	x2=BeamDataBase(i,4:6);
	for j=i:size(BeamDataBase,1)
		x3=BeamDataBase(j,1:3);
		x4=BeamDataBase(j,4:6);
		a=x2-x1;
		b=x4-x3;
		c=x3-x1;
		if i~=j 
			% are lines planar ?
			if dot((x1-x3),cross((x2-x1),(x4-x3)))==0
				q1=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
				q2=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
				v=0;
			else
				%position of minimum distance
				u1=a/norm(a);
				u2=b/norm(b);
				u3=cross(u1,u2)/norm(cross(u1,u2));
				e=dot((x1-x3),u3)*u3;
				v=norm(dot(c,cross(a,b)))/norm(cross(a,b));
				x3dash=x3+e;
				x4dash=x4+e;
				a=x2-x1;
				b=x4dash-x3dash;
				c=x3dash-x1;
				q1=x1+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
				e=dot((x3-x1),u3)*u3;
				x1dash=x1+e;
				x2dash=x2+e;				
				a=x2dash-x1dash;
				b=x4-x3;
				c=x3-x1dash;
				q2=x1dash+a*dot(cross(c,b),cross(a,b))/norm(cross(a,b))^2;
			end
			BalltoMinDist(1)=sqrt(sum((q1-x1).^2));
			BalltoMinDist(2)=sqrt(sum((q1-x2).^2));
			BalltoMinDist(3)=sqrt(sum((q1-x3).^2));
			BalltoMinDist(4)=sqrt(sum((q1-x4).^2));
			[a,b]=sort([x1(1) x2(1)]);
			[c,d]=sort([x3(1) x4(1)]);
			if q1(1)>a(1) & q1(1)<a(2) & q2(1)>c(1) & q2(1)<c(2)
				[a,b]=sort([x1(2) x2(2)]);
				[c,d]=sort([x3(2) x4(2)]);
				if q1(2)>a(1) & q1(2)<a(2) & q2(2)>c(1) & q2(2)<c(2)
					[a,b]=sort([x1(3) x2(3)]);
					[c,d]=sort([x3(3) x4(3)]);
					if q1(3)>a(1) & q1(3)<a(2) & q2(3)>c(1) & q2(3)<c(2)
						%%% Check if min distance is very close to joint and allow v>0.0053 if this is the case 
						if min(BalltoMinDist)<0.02
							if v < 0.0053
								V=V+0.0053-v;
							end
						else
							if v < Goal
								V=V+Goal-v;
							end
						end
					end
				end
			end			
		end
	end
end

Vbeams=V-Vvertical-Vjoints-Vjoints2;
Dist=V;


