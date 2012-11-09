% % SODES, symbolic differential equation solver
% William McGinnis
% 10/30/2012
%
% % Inputs
%   eqs:                {2nd order diffeq, 2nd order diffeq, ....}
%   free_var:           'x' or 'dx' or 'velocity', name of free variable.
%   t_init:             initial time, in seconds
%   t_final:            final time, in seconds
%   initial conditions: intial conditions, state vector, in columns
%
%   optional arguments:
%   'Integrator': options: 'ode45', 'ode23', 'ode15s', 'mexTrapz' (default ode45)
%   'TimeStep' : options: time step in seconds (default none)
%   'Plotting'  : options: ON or OFF (default OFF)
%   'Terminal'  : options: ON or OFF (default OFF)
%   'StateNames': options: Cell array of state names (default y1,y2,...,yn)
%   'StateUnits': options: Cell array of state units (default units)
%   'Options'   : send an odeset options object
% % Outputs
%   [t xs]
%   t: time vector
%   xs: output of all states in columns
% % Example
% Every feature:
% statenames={'theta 1', 'omega 1', 'theta 2', 'omega 2'};
% stateunits={'radians', 'radians/sec', 'radians', 'radians/sec'};
% initalconds=[rand;rand;rand;rand];
% [t xs]=sodes({ddq1s,ddq2s},'x',0,10,initialconds, ...
%      'Integrator', 'ode45','TimeStep', 0.1, 'Terminal', 'ON',...
%      'Plotting', 'ON','StateNames', statenames, 'StateUnits', stateunits);
%
% Bare Minimum:
% initalconds=[rand;rand;rand;rand];
% [t xs]=sodes({ddq1s,ddq2s},'x',0,10,initialconds);
%
%
%
% % Current Math Limitations:
%  Not all math will convert right at this time.  The main source of
%  problems centers around the differences in how C and MATLAB handle
%  exponentials.  C uses ^ as the bitwise XOR operator, MATLAB as pow.  To
%  get around this, there is some code to find any ^'s in the MATLAB code,
%  and back track to find the clause that it is raising to whatever power.
%  Right now the code will not convert correctly if:
%       -Anything is raised to a power that isnt -1<x<10 (one digit only)
%       -Any function (abs(x)^2, inv(x)^2, ect.) with the exception of:
%           -sin
%           -cos
%           -tan
%           -sec
%           -csc
%           -cot
%           -asin
%           -acos
%           -atan
%           -asec
%           -acsc
%           -acot
%           -sinh
%           -cosh
%           -tanh
%           -sech
%           -csch
%           -coth
%           -asinh
%           -acosh
%           -atanh
%           -asech
%           -acsch
%           -acoth
%           -exp
%           -log
%           -sqrt
%  Work will continue to account for more built in functions to be
%  correctly converted from MATLAB to C.
function [t, xs]=sodes(eqs, free_var, t_init, t_final, ic, varargin)
%% Input Handling/Housekeeping
filename=strcat(pwd,'\derivs.c');
fprintf('WORKING FILE: %s\n\n\n\n\n',filename);
num_states=size(eqs,2)*2;
Plotting=0;
Terminal=0;
TimeStep=0;
%1=ode45, 2=ode23, 3=ode15s
Integrator=1;
%% State Name Parsing
StateNames=cell(1,num_states);
for j=1:num_states
    StateNames{1,j}=sprintf('y%i',j);
end

for j=1:size(varargin,2)
    if ischar(varargin{1,j})
        if strcmp(varargin{1,j},'StateNames')
            if size(varargin{1,j+1},2)==num_states  || size(varargin{1,j+1},1)==num_states
                StateNames=varargin{1,j+1};
            else
                warning('Must have the same number of state names as states');
            end
        end
    end
end

%% State Units Parsing
StateUnits=cell(1,num_states);
for j=1:num_states
    StateUnits{1,j}='Units';
end

for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'StateUnits')
            if size(varargin{1,j+1},2)==num_states || size(varargin{1,j+1},1)==num_states
                StateUnits=varargin{1,j+1};
            else
                warning('Must have the same number of state units as states');
            end
        end
    end
end
%% Plotting Parsing
for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'Plotting')
            if strcmp(varargin{1,j+1},'ON')
                Plotting=1;
            else
                Plotting=0;
            end
        end
    end
end

%% Terminal Parsing
for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'Terminal')
            if strcmp(varargin{1,j+1},'ON')
                Terminal=1;
            else
                Terminal=0;
            end
        end
    end
end

%% Integrator Parsing
for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'Integrator')
            if strcmp(varargin{1,j+1},'ode45')
                Integrator=1;
            elseif strcmp(varargin{1,j+1},'ode23')
                Integrator=2;
            elseif strcmp(varargin{1,j+1},'ode15s')
                Integrator=3;
            else
                Integrator=4;
            end
        end
    end
end
%% TimeStep Parsing
for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'TimeStep')
            TimeStep=varargin{1,j+1};
        end
    end
end

%% Options Parsing
for j=1:size(varargin,2)
    if  ischar(varargin{1,j})
        if strcmp(varargin{1,j},'Options')
            ODEOptions=varargin{1,j+1};
        end
    end
end



%% Writing the MEX file for the integration itself.
if Integrator~=4
    fileid=fopen(filename, 'w');
    fprintf(fileid,'#include <math.h>\n');
    fprintf(fileid,'#include "mex.h"\n');
    fprintf(fileid,'#define	T_IN	prhs[0]\n');
    fprintf(fileid,'#define	Y_IN	prhs[1]\n');
    fprintf(fileid,'#define	YP_OUT	plhs[0]\n\n');
    fprintf(fileid,'#if !defined(MAX)\n');
    fprintf(fileid,'\t#define MAX(A, B)       ((A) > (B) ? (A) : (B))\n');
    fprintf(fileid,'#endif\n\n');
    fprintf(fileid,'#if !defined(MIN)\n');
    fprintf(fileid,'\t#define MIN(A, B)       ((A) < (B) ? (A) : (B))\n');
    fprintf(fileid,'#endif\n\n\n\n\n');
    fprintf(fileid,'static void yprime(double yp[],double *t,double y[])\n');
    fprintf(fileid,'{\n');
    
    %equations of motion
    state_count=1;
    for j=1:2:num_states
        fprintf(fileid,'\typ[%i]=',j-1);
        fprintf(fileid,'y[%i];\n',j);
        
        linetemp=char(eqs{1,state_count});
        if Terminal==1
            fprintf('*************Starting Point*****************\n')
            fprintf('%s\n',linetemp);
            fprintf('\n\n\n');
        end
        for k=1:num_states
            linetemp=strrep(linetemp,sprintf('%s(%i)',free_var,k),sprintf('y[%i]',k-1));
        end
        lineout=[];
        cursor=1;
        for k=1:length(linetemp)
            if strcmp('^', linetemp(k))
                carrot_pos=k;
                power_num=str2double(linetemp(k+1));
                if Terminal==1; fprintf('\n\nfound a power (x^%i) function, at position %i\n',power_num,carrot_pos); end;
                %Case: parenthetical clause
                if strcmp(linetemp(k-1),')')
                    if Terminal==1; fprintf('\tlooks like it is parenthetical...\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if strcmp(linetemp(m),'(') && open_pos==0
                            open_pos=m;
                            if Terminal==1; fprintf('\tfound the opening parentheses at %i\n',open_pos); end;
                        end
                    end
                    start_point=open_pos-1;
                    clause=linetemp(open_pos:k-1);
                    if Terminal==1; fprintf('\tthe base clause is: %s\n',clause); end;
                    if open_pos~=0 && open_pos>=4 && open_pos~=' '
                        if Terminal==1; fprintf('\tchecking for a trig function or something...\n'); end;
                        pre_clause=linetemp(open_pos-3:open_pos-1);
                        if strcmp(pre_clause,'sin') || strcmp(pre_clause,'cos') ...
                                || strcmp(pre_clause,'tan') || strcmp(pre_clause,'sec')...
                                || strcmp(pre_clause,'csc') || strcmp(pre_clause,'cot')...
                                || strcmp(pre_clause,'exp') || strcmp(pre_clause,'log')
                            if Terminal==1; fprintf('\tfound one!\n'); end;
                            if strcmp(open_pos-4,'a')
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            else
                                clause=strcat(pre_clause,clause);
                                start_point=open_pos-4;
                            end
                            if Terminal==1; fprintf('\tnew clause: %s\n',clause); end;
                        elseif strcmp(linetemp(open_pos-4:open_pos-1),'sinh') || strcmp(linetemp(open_pos-4:open_pos-1),'cosh')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'tanh') || strcmp(linetemp(open_pos-4:open_pos-1),'sech')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'csch') || strcmp(linetemp(open_pos-4:open_pos-1),'coth')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'sqrt') ||  strcmp(linetemp(open_pos-4:open_pos-1),'ceil')
                            if strcmp(open_pos-5,'a')
                                clause=strcat(linetemp(open_pos-5:open_pos-1),clause);
                                start_point=open_pos-6;
                            else
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            end
                        end
                    end
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:start_point));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                else
                    if Terminal==1; fprintf('\tthis isn''t parenthetical\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if(strcmp(linetemp(m), ' ') || strcmp(linetemp(m), '*') ...
                                || strcmp(linetemp(m), '+') || strcmp(linetemp(m), '-') ...
                                || strcmp(linetemp(m), '/'))
                            if open_pos==0
                                open_pos=m;
                                if Terminal==1; fprintf('\tBeginning of clause found at %i\n', m); end;
                            end
                        end
                    end
                    clause=linetemp(open_pos+1:k-1);
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:open_pos));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                end
            end
        end
        lineout=strcat(lineout, linetemp(cursor:end));
        fprintf(fileid,'\typ[%i]=%s;\n',j,lineout);
        state_count=state_count+1;
    end
    
    
    
    fprintf(fileid,'}\n\n\n\n');
    fprintf(fileid,'void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )\n');
    fprintf(fileid,'{\n');
    fprintf(fileid,'\tdouble *yp;\n');
    fprintf(fileid,'\tdouble *t,*y;\n');
    fprintf(fileid,'\tunsigned int m,n;\n\n');
    fprintf(fileid,'\tif (nrhs != 2) {\n');
    fprintf(fileid,'\t\tmexErrMsgTxt("Two input arguments required.");\n');
    fprintf(fileid,'\t} else if (nlhs > 1) {\n');
    fprintf(fileid,'\tmexErrMsgTxt("Too many output arguments.");\n');
    fprintf(fileid,'\t}\n');
    fprintf(fileid,'\t/* Check the dimensions of Y.  Y can be %i X 1 or 1 X %i. */ \n',num_states, num_states);
    fprintf(fileid,'\tm = mxGetM(Y_IN);\n');
    fprintf(fileid,'\tn = mxGetN(Y_IN);\n');
    fprintf(fileid,'\tif (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) ||\n');
    fprintf(fileid,'\t\t(MAX(m,n) != %i) || (MIN(m,n) != 1)) {\n',num_states);
    fprintf(fileid,'\t\tmexErrMsgTxt("mexODEFun requires that Y be a %i x 1vector.");\n',num_states);
    fprintf(fileid,'\t}\n');
    fprintf(fileid,'\t/* Create a matrix for the return argument */ \n');
    fprintf(fileid,'\tYP_OUT = mxCreateDoubleMatrix(m, n, mxREAL);\n');
    fprintf(fileid,'\t/* Assign pointers to the various parameters */ \n');
    fprintf(fileid,'\typ = mxGetPr(YP_OUT);\n');
    fprintf(fileid,'\tt = mxGetPr(T_IN);\n');
    fprintf(fileid,'\ty = mxGetPr(Y_IN);\n');
    fprintf(fileid,'\t/* Do the actual computations in a subroutine */ \n');
    fprintf(fileid,'\typrime(yp,t,y);\n');
    fprintf(fileid,'\t/* Return YP_OUT to MATLAB environment */ \n');
    fprintf(fileid,'\treturn;\n');
    fprintf(fileid,'}\n');
    fclose(fileid);
    
    %% Compiling the MEX file
    mex derivs.c
    
    
    %% The integration itself
    if TimeStep==0
        tspan=[t_init t_final];
    else
        tspan=t_init:TimeStep:t_final;
    end
    y0=ic;
    if Integrator==1
        if exist('ODEOptions','var')~=0
            [t,xs]= ode45('derivs', tspan, y0, ODEOptions);
        else
            [t,xs]= ode45('derivs', tspan, y0);
        end
    elseif Integrator==2
        if exist('ODEOptions','var')~=0
            [t,xs]= ode23('derivs', tspan, y0, ODEOptions);
        else
            [t,xs]= ode23('derivs', tspan, y0);
        end
    else
        if exist('ODEOptions','var')~=0
            [t,xs]= ode15s('derivs', tspan, y0, ODEOptions);
        else
            [t,xs]= ode15s('derivs', tspan, y0);
        end
    end
else
    fileid=fopen(filename, 'w');
    fprintf(fileid,'#include <math.h>\n');
    fprintf(fileid,'#include "mex.h"\n');
    fprintf(fileid,'#define	T_IN	prhs[0]\n');
    fprintf(fileid,'#define	Y_IN	prhs[1]\n');
    fprintf(fileid,'#define	YP_OUT	plhs[0]\n\n');
    fprintf(fileid,'#if !defined(MAX)\n');
    fprintf(fileid,'\t#define MAX(A, B)       ((A) > (B) ? (A) : (B))\n');
    fprintf(fileid,'#endif\n\n');
    fprintf(fileid,'#if !defined(MIN)\n');
    fprintf(fileid,'\t#define MIN(A, B)       ((A) < (B) ? (A) : (B))\n');
    fprintf(fileid,'#endif\n\n\n\n\n');
    fprintf(fileid,'static void yprime(double yp[],double *t,double y[])\n');
    fprintf(fileid,'{\n');
    fprintf(fileid,'\tdouble k1_0, k1_1, k1_2, k1_3,k2_0, k2_1, k2_2, k2_3, y1_1, y1_2, y1_0, y1_3;\n');
    fprintf(fileid,'\tdouble yp_k[4]={0,0,0,0};\n');
    fprintf(fileid,'\tdouble x0=%9.5f;\n',t_init);
    fprintf(fileid,'\tdouble tf=%9.5f;\n',t_final);
    if TimeStep==0
        fprintf(fileid,'\tdouble h=0.1;\n');
    else
        fprintf(fileid,'\tdouble h=%9.5f;\n', TimeStep);
    end
    fprintf(fileid,'\tint io;\n');
    fprintf(fileid,'\tFILE *file;\n');
    fprintf(fileid,'\tfile=fopen("integration_data.dat","w");\n');
    fprintf(fileid,'\tio=fprintf(file, "%%f\\t%%f\\t%%f\\t%%f\\n",y[0],y[1],y[2],y[3]);\n');
    
    
    
    fprintf(fileid,'\tfor(; x0<tf; )\n');
    fprintf(fileid,'\t{\n');
    
    fprintf(fileid,'\t\tx0=x0+h;\n');
    
    %equations of motion
    state_count=1;
    for j=1:2:num_states
        fprintf(fileid,'\t\typ[%i]=',j-1);
        fprintf(fileid,'y[%i];\n',j);
        
        linetemp=char(eqs{1,state_count});
        if Terminal==1
            fprintf('*************Starting Point*****************\n')
            fprintf('%s\n',linetemp);
            fprintf('\n\n\n');
        end
        for k=1:num_states
            linetemp=strrep(linetemp,sprintf('%s(%i)',free_var,k),sprintf('y[%i]',k-1));
        end
        lineout=[];
        cursor=1;
        for k=1:length(linetemp)
            if strcmp('^', linetemp(k))
                carrot_pos=k;
                power_num=str2double(linetemp(k+1));
                if Terminal==1; fprintf('\n\nfound a power (x^%i) function, at position %i\n',power_num,carrot_pos); end;
                %Case: parenthetical clause
                if strcmp(linetemp(k-1),')')
                    if Terminal==1; fprintf('\tlooks like it is parenthetical...\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if strcmp(linetemp(m),'(') && open_pos==0
                            open_pos=m;
                            if Terminal==1; fprintf('\tfound the opening parentheses at %i\n',open_pos); end;
                        end
                    end
                    start_point=open_pos-1;
                    clause=linetemp(open_pos:k-1);
                    if Terminal==1; fprintf('\tthe base clause is: %s\n',clause); end;
                    if open_pos~=0 && open_pos>=4 && open_pos~=' '
                        if Terminal==1; fprintf('\tchecking for a trig function or something...\n'); end;
                        pre_clause=linetemp(open_pos-3:open_pos-1);
                        if strcmp(pre_clause,'sin') || strcmp(pre_clause,'cos') ...
                                || strcmp(pre_clause,'tan') || strcmp(pre_clause,'sec')...
                                || strcmp(pre_clause,'csc') || strcmp(pre_clause,'cot')...
                                || strcmp(pre_clause,'exp') || strcmp(pre_clause,'log')
                            if Terminal==1; fprintf('\tfound one!\n'); end;
                            if strcmp(open_pos-4,'a')
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            else
                                clause=strcat(pre_clause,clause);
                                start_point=open_pos-4;
                            end
                            if Terminal==1; fprintf('\tnew clause: %s\n',clause); end;
                        elseif strcmp(linetemp(open_pos-4:open_pos-1),'sinh') || strcmp(linetemp(open_pos-4:open_pos-1),'cosh')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'tanh') || strcmp(linetemp(open_pos-4:open_pos-1),'sech')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'csch') || strcmp(linetemp(open_pos-4:open_pos-1),'coth')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'sqrt') ||  strcmp(linetemp(open_pos-4:open_pos-1),'ceil')
                            if strcmp(open_pos-5,'a')
                                clause=strcat(linetemp(open_pos-5:open_pos-1),clause);
                                start_point=open_pos-6;
                            else
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            end
                        end
                    end
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:start_point));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                else
                    if Terminal==1; fprintf('\tthis isn''t parenthetical\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if(strcmp(linetemp(m), ' ') || strcmp(linetemp(m), '*') ...
                                || strcmp(linetemp(m), '+') || strcmp(linetemp(m), '-') ...
                                || strcmp(linetemp(m), '/'))
                            if open_pos==0
                                open_pos=m;
                                if Terminal==1; fprintf('\tBeginning of clause found at %i\n', m); end;
                            end
                        end
                    end
                    clause=linetemp(open_pos+1:k-1);
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:open_pos));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                end
            end
        end
        lineout=strcat(lineout, linetemp(cursor:end));
        fprintf(fileid,'\t\typ[%i]=%s;\n',j,lineout);
        state_count=state_count+1;
    end
    
    fprintf(fileid,'\t\tk1_0 = h * (yp[0]);\n');
    fprintf(fileid,'\t\tk1_1 = h * (yp[1]);\n');
    fprintf(fileid,'\t\tk1_2 = h * (yp[2]);\n');
    fprintf(fileid,'\t\tk1_3 = h * (yp[3]);\n');
    
    fprintf(fileid,'\t\typ_k[0]=(yp[0]) + k1_0;\n');
    fprintf(fileid,'\t\typ_k[1]=(yp[1]) + k1_1;\n');
    fprintf(fileid,'\t\typ_k[2]=(yp[2]) + k1_2;\n');
    fprintf(fileid,'\t\typ_k[3]=(yp[3]) + k1_3;\n');
    
    fprintf(fileid,'\t\tx0=x0+h;\n');
    
    
    %equations of motion
    state_count=1;
    for j=1:2:num_states
        fprintf(fileid,'\t\typ[%i]=',j-1);
        fprintf(fileid,'y[%i];\n',j);
        
        linetemp=char(eqs{1,state_count});
        if Terminal==1
            fprintf('*************Starting Point*****************\n')
            fprintf('%s\n',linetemp);
            fprintf('\n\n\n');
        end
        for k=1:num_states
            linetemp=strrep(linetemp,sprintf('%s(%i)',free_var,k),sprintf('y[%i]',k-1));
        end
        lineout=[];
        cursor=1;
        for k=1:length(linetemp)
            if strcmp('^', linetemp(k))
                carrot_pos=k;
                power_num=str2double(linetemp(k+1));
                if Terminal==1; fprintf('\n\nfound a power (x^%i) function, at position %i\n',power_num,carrot_pos); end;
                %Case: parenthetical clause
                if strcmp(linetemp(k-1),')')
                    if Terminal==1; fprintf('\tlooks like it is parenthetical...\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if strcmp(linetemp(m),'(') && open_pos==0
                            open_pos=m;
                            if Terminal==1; fprintf('\tfound the opening parentheses at %i\n',open_pos); end;
                        end
                    end
                    start_point=open_pos-1;
                    clause=linetemp(open_pos:k-1);
                    if Terminal==1; fprintf('\tthe base clause is: %s\n',clause); end;
                    if open_pos~=0 && open_pos>=4 && open_pos~=' '
                        if Terminal==1; fprintf('\tchecking for a trig function or something...\n'); end;
                        pre_clause=linetemp(open_pos-3:open_pos-1);
                        if strcmp(pre_clause,'sin') || strcmp(pre_clause,'cos') ...
                                || strcmp(pre_clause,'tan') || strcmp(pre_clause,'sec')...
                                || strcmp(pre_clause,'csc') || strcmp(pre_clause,'cot')...
                                || strcmp(pre_clause,'exp') || strcmp(pre_clause,'log')
                            if Terminal==1; fprintf('\tfound one!\n'); end;
                            if strcmp(open_pos-4,'a')
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            else
                                clause=strcat(pre_clause,clause);
                                start_point=open_pos-4;
                            end
                            if Terminal==1; fprintf('\tnew clause: %s\n',clause); end;
                        elseif strcmp(linetemp(open_pos-4:open_pos-1),'sinh') || strcmp(linetemp(open_pos-4:open_pos-1),'cosh')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'tanh') || strcmp(linetemp(open_pos-4:open_pos-1),'sech')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'csch') || strcmp(linetemp(open_pos-4:open_pos-1),'coth')...
                                || strcmp(linetemp(open_pos-4:open_pos-1),'sqrt') ||  strcmp(linetemp(open_pos-4:open_pos-1),'ceil')
                            if strcmp(open_pos-5,'a')
                                clause=strcat(linetemp(open_pos-5:open_pos-1),clause);
                                start_point=open_pos-6;
                            else
                                clause=strcat(linetemp(open_pos-4:open_pos-1),clause);
                                start_point=open_pos-5;
                            end
                        end
                    end
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:start_point));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                else
                    if Terminal==1; fprintf('\tthis isn''t parenthetical\n'); end;
                    open_pos=0;
                    for m=k-1:-1:1
                        if(strcmp(linetemp(m), ' ') || strcmp(linetemp(m), '*') ...
                                || strcmp(linetemp(m), '+') || strcmp(linetemp(m), '-') ...
                                || strcmp(linetemp(m), '/'))
                            if open_pos==0
                                open_pos=m;
                                if Terminal==1; fprintf('\tBeginning of clause found at %i\n', m); end;
                            end
                        end
                    end
                    clause=linetemp(open_pos+1:k-1);
                    if Terminal==1; fprintf('\t\tExpanding it out....\n'); end;
                    lineout=strcat(lineout,linetemp(cursor:open_pos));
                    expclause=clause;
                    for m=power_num-1
                        expclause=strcat(expclause,strcat('*',clause));
                    end
                    lineout=strcat(lineout,expclause);
                    if Terminal==1; fprintf('\t\texpanded clause: %s\n',expclause); end;
                    if Terminal==1; fprintf('\t\tEquation so far: %s\n',lineout); end;
                    cursor=carrot_pos+2;
                end
            end
        end
        lineout=strcat(lineout, linetemp(cursor:end));
        fprintf(fileid,'\t\typ[%i]=%s;\n',j,lineout);
        state_count=state_count+1;
    end
    
    fprintf(fileid,'\t\tx0=x0-h;\n');
    
    fprintf(fileid,'\t\tk2_0 = h * (yp[0]);\n');
    fprintf(fileid,'\t\tk2_1 = h * (yp[1]);\n');
    fprintf(fileid,'\t\tk2_2 = h * (yp[2]);\n');
    fprintf(fileid,'\t\tk2_3 = h * (yp[3]);\n');
    
    fprintf(fileid,'\t\ty1_0 = (y[0]) + ( k1_0 + k2_0)/2;\n');
    fprintf(fileid,'\t\ty1_1 = (y[1]) + ( k1_1 + k2_1)/2;\n');
    fprintf(fileid,'\t\ty1_2 = (y[2]) + ( k1_2 + k2_2)/2;\n');
    fprintf(fileid,'\t\ty1_3 = (y[3]) + ( k1_3 + k2_3)/2;\n');
    
    fprintf(fileid,'\t\ty[0]=y1_0;\n');
    fprintf(fileid,'\t\ty[1]=y1_1;\n');
    fprintf(fileid,'\t\ty[2]=y1_2;\n');
    fprintf(fileid,'\t\ty[3]=y1_3;\n');
    fprintf(fileid,'\t\tio=fprintf(file, "%%f\\t%%f\\t%%f\\t%%f\\n",y[0],y[1],y[2],y[3]);\n');
	fprintf(fileid,'\t\t}\n');
    fprintf(fileid,'\tio=fclose(file);\n');
    fprintf(fileid,'\typ[0]=y[0];\n');
    fprintf(fileid,'\typ[1]=y[1];\n');
    fprintf(fileid,'\typ[2]=y[2];\n');
    fprintf(fileid,'\typ[3]=y[3];\n');
    
    
    
    
    fprintf(fileid,'}\n\n\n\n');
    fprintf(fileid,'void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )\n');
    fprintf(fileid,'{\n');
    fprintf(fileid,'\tdouble *yp;\n');
    fprintf(fileid,'\tdouble *t,*y;\n');
    fprintf(fileid,'\tunsigned int m,n;\n\n');
    fprintf(fileid,'\tif (nrhs != 2) {\n');
    fprintf(fileid,'\t\tmexErrMsgTxt("Two input arguments required.");\n');
    fprintf(fileid,'\t} else if (nlhs > 1) {\n');
    fprintf(fileid,'\tmexErrMsgTxt("Too many output arguments.");\n');
    fprintf(fileid,'\t}\n');
    fprintf(fileid,'\t/* Check the dimensions of Y.  Y can be %i X 1 or 1 X %i. */ \n',num_states, num_states);
    fprintf(fileid,'\tm = mxGetM(Y_IN);\n');
    fprintf(fileid,'\tn = mxGetN(Y_IN);\n');
    fprintf(fileid,'\tif (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) ||\n');
    fprintf(fileid,'\t\t(MAX(m,n) != %i) || (MIN(m,n) != 1)) {\n',num_states);
    fprintf(fileid,'\t\tmexErrMsgTxt("mexODEFun requires that Y be a %i x 1vector.");\n',num_states);
    fprintf(fileid,'\t}\n');
    fprintf(fileid,'\t/* Create a matrix for the return argument */ \n');
    fprintf(fileid,'\tYP_OUT = mxCreateDoubleMatrix(m, n, mxREAL);\n');
    fprintf(fileid,'\t/* Assign pointers to the various parameters */ \n');
    fprintf(fileid,'\typ = mxGetPr(YP_OUT);\n');
    fprintf(fileid,'\tt = mxGetPr(T_IN);\n');
    fprintf(fileid,'\ty = mxGetPr(Y_IN);\n');
    fprintf(fileid,'\t/* Do the actual computations in a subroutine */ \n');
    fprintf(fileid,'\typrime(yp,t,y);\n');
    fprintf(fileid,'\t/* Return YP_OUT to MATLAB environment */ \n');
    fprintf(fileid,'\treturn;\n');
    fprintf(fileid,'}\n');
    fclose(fileid);
    
    mex derivs.c
    
    temp=derivs(t_init, ic);
        
    xs=load('integration_data.dat');
    t=linspace(t_init,t_final, size(xs,1));
end


if Plotting==1
    figure;
    for j=1:num_states
        subplot(num_states,1,j)
        plot(t,xs(:,j),'b')
        xlabel('Time, (Seconds)');
        ylabel(StateUnits{1,j});
        title(StateNames{1,j});
    end
end