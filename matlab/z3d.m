classdef z3d<handle
% Y=z3d(H,X) where H \in R^{3m} and X=[x1min x1max;x2min x2max;...]
% constructs 3D zonotope Y from matrix of generators H with m columns and 
% matrix X of interval constraints on each coordinate 
% m<=64 as we use bitwise operations
% The zonotope initially is constructed in the form of linear inequalities
%
% Y.extreme_points() calculates its vertices
% Y.draw(color) displays it (color is the same char argument as in 'plot')
%               Note that you cannot draw zonotope without calculating 
%               its extreme points
% Y.facets_num() returns number of facets (actually, only half of them as 
%                the zonotope is zero-symmetric)
% [d,c]=Y.facet(i) returns linear inequality for i-th facet so that
%       <d,x><=c
% y=Y.vertices(i) returns vertices as columns of y for the corresponding 
%                 i-th facet
%
% Author: Max Demenkov, Institute of Control Sciences, Moscow, 2013-2016
%
    properties
     N=0;
     F=[];h=[];
     Vy=[];Vu=[];
     y0=[];u0=[];
     f2v=[];list=[];
     B0=[];umax=[];
     tol=1000*eps;% Tweak it in case of problems
    end
    
    methods(Static)
        function [flag,B0,U0,y0,u0]=get_symmetric(B,U)
        % Makes the zonotope zero-centered
                    [m,n]=size(B);flag=1;
                    B0=[];U0=[];y0=zeros(m,1);
                    for i=1:n
                        if U(i,1)>U(i,2)
                           return;
                        elseif U(i,1)==U(i,2)
                           y0=y0+B(:,i)*U(i,1);
                        elseif z3d.strong_zero(B(:,i),1e-4)
                        else
                           u0=(U(i,1)+U(i,2))/2;
                           U0=[U0;U(i,1)-u0 U(i,2)-u0];
                           B0=[B0,B(:,i)];y0=y0+B(:,i)*u0;
                        end
                    end
                    flag=0;
        end

                function N=face_num(controls)

                N=0;
                for i=1:(controls-1)
                    N=N+i;
                end
                end
                
                function out=arrange(list,Y,tol)

                in=list;SV=Y(:,list);
                out=in(1);cur=out;in(1)=[];
                while ~isempty(in)
                      n1=length(in);
                      for i=1:length(in)
                          d=z3d.norm_vec(cross(Y(:,cur),Y(:,in(i))))';
                          if all(d*SV<=tol) || all(d*SV>=-tol)
                             cur=in(i);out=[out cur];in(i)=[];
                             break;
                          end
                      end
                      n2=length(in);
                      if n1==n2, error('Arrange failure'); end 
                end
                end
                
                function f2v=real_vertices(VE,n,list)

                fn=length(list);k=length(VE);
                CM=zeros(fn,k);
                for i=1:fn
                    for j=1:length(list{i})
                       num=list{i}(j);vert=VE(num);
                       CM(i,num)=1;
                       ind=find(VE==bitxor(2^n-1,vert));
                       CM(i,ind)=-1;
                    end  
                end
                remove=find(sum(abs(CM))<3);
                for i=remove
                    CM(:,i)=zeros(fn,1);
                end
                for i=1:fn
                    f2v{i}=find(CM(i,:)>0);
                end
                end
                
                function [VE,f2v]=unique_vertices(list)

                fn=length(list);
                VE(1,1)=list{1}(1);k=1;f2v{1}(1)=1;
                for i=1:fn
                    for j=1:length(list{i})
                        vert=list{i}(j);
                        ind=find(VE==vert);
                        if isempty(ind)
                           k=k+1;VE(k,1)=vert;f2v{i}(j)=k;
                        else
                           f2v{i}(j)=ind;
                        end
                    end
                end
                end
                
                function [keep,remove]=unique_faces(F1,tol)

                [m,n]=size(F1);keep=1;remove=[];
                F2=F1(1,:);k=1;
                for i=2:m
                    add_this=1;
                    for j=1:k
                        if z3d.zero(cross(F2(j,:),F1(i,:)),tol)
                           add_this=0;break;
                        end  
                    end
                    if add_this
                       k=k+1;F2(k,:)=F1(i,:);keep=[keep;i];
                    else
                       remove=[remove;i];
                    end
                end
                end
                
                function E=code(d,tol)

                ind=find(abs(d)<=tol);
                s=length(ind);
                u=0;E=[];
                for i=find(d>tol)
                    u=bitset(u,i,1);
                end
                for i=0:(2^s-1)
                    for j=1:s
                        u=bitset(u,ind(j),bitget(i,j));
                    end
                    E=[E,u];
                end
                end
                
                function Y=decode_Y(E,B)

                [m,n]=size(B);k=length(E);Y=zeros(m,k);
                for j=1:k
                    y=zeros(m,1);
                    for i=1:n
                        if bitget(E(j),i),dy=B(:,i);else,dy=-B(:,i);end
                        y=y+dy;
                    end
                    Y(:,j)=y;
                end
                end
                
                function U=decode_U(E,umax)

                n=length(umax);m=length(E);U=zeros(n,m);
                for j=1:m
                    v=zeros(n,1);
                    for i=1:n
                        if bitget(E(j),i),v(i)=umax(i);else,v(i)=-umax(i);end
                    end
                    U(:,j)=v;
                end
                end
                
                function S2=norm_vec(S1)

                [m,n]=size(S1);S2=zeros(m,n);
                for i=1:n
                    d=S1(:,i);
                    S2(:,i)=d./sqrt(d'*d);
                end
                end
                
                function B2=b_real(B1,umax)

                [m,n]=size(B1);B2=zeros(m,n);
                for i=1:n
                    B2(:,i)=B1(:,i)*umax(i);
                end
                end
                
                function flag=zero(vec,tol)

                flag=max(abs(vec))<=tol;
                end

                function flag=strong_zero(vec,tol)

                flag=max(abs(vec))<tol;
                end
                
                function [S,rhs]=normalize(S,rhs)

                [ms,ns]=size(S);
                for k=1:ms
                    if (abs(rhs(k))>eps)
                       S(k,:)=S(k,:)./abs(rhs(k));
                       rhs(k)=sign(rhs(k));
                    end
                end
                end
                
                function S2=addcol(S1,v)
                [m,n]=size(S1);
                for i=1:n
                    S2(:,i)=S1(:,i)+v;
                end
                end
                
                
    end
    
    methods
        
        function self=z3d(B,U,varargin)
        % Constructs zonotope as a system of linear inequalities
            if nargin>2,tol=varargin{3}; else tol=self.tol; end
            [flag,B0,U0,self.y0,self.u0]=z3d.get_symmetric(B,U);
            if flag, error('Wrong box limits'); end
            [m,n]=size(B0);umax=U0(:,2);self.umax=umax;
            fn=z3d.face_num(n);
                BN=z3d.norm_vec(B0);
                BR=z3d.b_real(B0,umax);
                ind=0;F=zeros(fn,m);h=zeros(fn,1);list=cell(fn,1);
                for i=1:n
                    for j=(i+1):n
                        face=cross(BN(:,i),BN(:,j))';
                        if ~z3d.zero(face,tol)
                           face=z3d.norm_vec(face')';
                           d=z3d.norm_vec((face*BN)')';
                           E=z3d.code(d,tol);y=z3d.decode_Y(E(1),BR);
                           val=face*y;
                           ind=ind+1;
                           F(ind,:)=face;h(ind,1)=val;list{ind}=E;
                        end
                    end  
                end
                F=F(1:ind,:);h=h(1:ind);list=list(1:ind);
                [keep,remove]=z3d.unique_faces(F,tol);
                fn=length(keep);
                F=F(keep,:);h=h(keep);list(remove)=[];
                [self.F,self.h]=z3d.normalize(F,h);
                self.N=2*fn;self.list=list;
                self.B0=B0;
              
        end
        
        function extreme_points(self,varargin)
        % Computes extreme points
                if nargin>1,tol=varargin{2}; else tol=self.tol; end
                BR=z3d.b_real(self.B0,self.umax);
                [VE,f2v]=z3d.unique_vertices(self.list);
                n=size(self.B0,2);
                f2v=z3d.real_vertices(VE,n,f2v);
                Y=z3d.decode_Y(VE,BR);U=z3d.decode_U(VE,self.umax);
                YN=z3d.norm_vec(Y);
                for i=1:length(f2v)
                    f2v{i}=z3d.arrange(f2v{i},YN,self.tol);
                end
                self.Vy=Y;self.Vu=U;
                self.f2v=f2v;
        end         
             
        function N=facets_num(self)
            N=self.N;
        end
        
        function [d,c]=facet(self,num)
            d=self.F(num,:);c=self.h(num)+d*self.y0;
        end
        
        function y=vertices(self,num)
            ind=self.f2v{num};
            y=self.Vy(:,ind);
            y=z3d.addcol(y,self.y0);
        end
        
        function draw(self,col)
        % Draws zonotope in 3D
               figure(); 
               hold on;material metal;
               f2v=self.f2v;
               Y1=z3d.addcol(self.Vy,self.y0);
               Y2=z3d.addcol(-self.Vy,self.y0);
               f2vm=zeros(length(f2v),length(f2v{1}));
               for i=1:length(f2v)
                   f2vm(i,:)=f2v{i};
               end
               Y=[Y1';Y2'];
               f2vm=[f2vm;ones(size(f2vm))*size(Y1,2)+f2vm];
               h=patch('vertices',Y,'faces',f2vm,'facecolor',col,...
                       'edgecolor','k','facealpha',1);
               facecolor = repmat([1 1 1],length(f2vm),1);
               set(h,'FaceVertexCData',facecolor);
               alpha(0.8);lighting phong;
               view([20 15]);
               camlight headlight;grid on;
               axis tight;
               
              
        end
        

    end
    
end
