% Weighted linear fit by chi^2 minimization.
%
% Last update: 12.07.13
%
% The syntax is
%   [parameters errors resultfunction Kovmatrix chi2/dof prob]=wlinfit(x,y,yerr,funccell,options)
%
%
% x,y,yerr are the data variables as vectors (arbitrary but equal orientation).
%
% funccell is a cell containing the functions which are used for linear
% fitting. For Example funccell = {@(x) sin(x), @(x) sin(2*x), @(x) sin(3*x)}
% For polynomial fit one can either add polynomial functions to the the
% funccell or specify the degree of the polynomial with the option polyfit
% (wlinfit(...,'polyfit',2) for second order).
%
% Important: All text input is in latex
%
% Options:
%
% Options are entered like: wlinfit(x,y,yerr,funccell,'optionname',value,'optionname2',value2)
%
%
%               
% polyfit: natural number                       Uses a polynomial of order
%                                               of the entered number to
%                                               fit. Additional functions
%                                               can be entered in the
%                                               funccell. If only a
%                                               polynomaial should be used
%                                               enter empty funccell ({}).
%                                               The functions in the
%                                               funccell will be listed
%                                               first. Default: -1 (no polynomial)
%                                                       
% 
% label:{'xaxis' 'yaxis' 'varname1'...}         Defines label of axis and arguments
%                                               Default: {'xaxis' 'yaxis' '$c_1$' '$c_2$'}
%
%
% plot: 'on' , true or 'off' , false            Defines whether plot is created. Plotting takes a lot of time)
%                                               Default: true
%
%
% header:   'headerstring' or                   Defines header in the legend. Use str if only one coloumn or
%           {'first col' 'second col'}          cell if multiple coloumns. Use '' for empty coloumn.
%                                               Default: 'Nonlinear Fit'
%
% errprec                                       Defines how many numerics
%                                               of value and error are
%                                               printed in the legend.
%                                               Default: 2
%
%
% axis: [xmin xmax ymin ymax]			Define axis for plot
%						Default: calculated from data
%
% format: [ xsize(cm) ysize(cm) ]               Defines size of the plots.
%                                               Default: [16 14]
%
%
% position: [xpos ypos]                         Defines position of legend. Use stupid input to swap legend.
%                                               Default: northwest or northeast depending on data
%
%
% grid: 'off' or 'y' or 'x' or 'on'             Defines whether and which kind of grid is plotted
%                                               Default: false
%
%
% printchi: 'off' or 'on'                       Defines whether chi^2/dof
%                                               is printed in legend
%                                               Default: true
%
%
% print: 'on' , true , 1 or 'off' , false , 0   Defines whether information is printed on command window
%                                               Default: true
%
% normalitytests:                               Defines whether
% 'on' , true , 1 or 'off' , false , 0          normalitytests are printed
%                                               on the screen
% 
% contourplt: matrices of 1x2 numbers			For every matrix a correlation
%												contourplot between the defined variables is created 
%												Default: no contourplot
%												Example: {[1 2],[1 3]}
%												creates correlation plots
%												of variable (1 and 2) and (1
%												and 3)
%
% Written by: Franziska Flegel & Jannick Weisshaupt (jannick@physik.hu-berlin.de)
% Paramtext by: Franziska Flegel & Jannick Weisshaupt
%
% Feel free to share, change or do whatever you want with this script.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%




function [beta,betaerr,funcout,Kov,chi2dof,prob] = wlinfit(x,y,yerr,funccell,varargin)


%%%%%%%%%%%%%%%%%%%%%%% Testing %%%%%%%%%%%%%%%%%%%


if any(~isnumeric(x)) || any(~isnumeric(y)) || any(~isnumeric(yerr))
    error('There are non numerical entries in your data (x,y or yerr)')
end

if ~iscell(funccell)
    error('Functions were not inserted as cell. If you use npolyfit only enter empty cell ({}).If you use only one function put it into {}')
end

if length(yerr)==1
    yerr=yerr *ones(size(y));
end

if any(size(x)~=size(y) )|| any(size(x)~=size(yerr) )
    error('Input Data must have same size')
end

if any(yerr==0)
    error('Errors with value 0 are not allowed.')
end


if ~any(y~=y(1))
   error('All y-values are equal or only one y value entered') 
end

if ~any(x~=x(1))
   error('All x-values are equal or only one x value entered') 
end

if any(y+yerr==y)
    warning('MATLAB:wnonlinfit:small_error_warning','Errors are smaller than machine accuracy. This might lead to very strange results')
end

[nx1 nx2] = size(x);
if nx2 ==1 && nx1>1
    x=x';
    y=y';
    yerr=yerr';
end
% Optional arguments

p = inputParser;

p.addParamValue('plot',true);
p.addParamValue('label',[]);
p.addParamValue('position',[])
p.addParamValue('format',[16 14])
p.addParamValue('grid',false)
p.addParamValue('print',false)
p.addParamValue('header',{'Linear Fit'})
p.addParamValue('printchi',true)
p.addParamValue('axis',[]);
p.addParamValue('errprec',2);
p.addParamValue('rescalex',1,@(x) isnumeric(x) && length(x)==1);
p.addParamValue('rescaley',1,@(x) isnumeric(x) && length(x)==1);
p.addParamValue('plotmode',[]);
p.addParamValue('polyfit',-1);
p.addParamValue('font',[]);
p.addParamValue('contourplot',{});
p.addParamValue('normalitytests',false);

p.parse(varargin{:});

stringcell=p.Results.label;
position=p.Results.position;
plotbool=p.Results.plot;
gridbool=p.Results.grid;
format=p.Results.format;
textheader=p.Results.header;
printbool=p.Results.print;
chibool=p.Results.printchi;
axisin=p.Results.axis;
errprec=p.Results.errprec;
font=p.Results.font;
rescalex=p.Results.rescalex;
rescaley=p.Results.rescaley;
plotmodestr=p.Results.plotmode;
npolyfit=p.Results.polyfit+1;
contourcell=p.Results.contourplot;
normalitytestbool=p.Results.normalitytests;

%%%%%% Convert optional input
if islogical(gridbool)
    if gridbool
        xgridstr='on';
        ygridstr='on';
    else
        xgridstr='off';
        ygridstr='off';
    end
end
if ischar(gridbool)
    if any(strcmpi(gridbool,{'y','yaxis','y-axis','yachse'}))
        xgridstr='off';
        ygridstr='on';
    elseif any(strcmpi(gridbool,{'x','xaxis','x-axis','xachse'}))
        xgridstr='on';
        ygridstr='off';
    elseif any(strcmpi(gridbool,{'off' 'n' 'no' 'non' 'neither'}))
        xgridstr='off';
        ygridstr='off';
    elseif any(strcmpi(gridbool,{'y','yes','ja','oui','si','on','true'}))
        xgridstr='on';
        ygridstr='on';
    else
        warning('MATLAB:locicalvar:inaprinput','Inapropriate input for grid. Standard value off is used')
        xgridstr='off';
        ygridstr='off';
    end
end

N=length(funccell)+npolyfit;
plotbool=logicalvar(plotbool,'plot');
chibool=logicalvar(chibool,'printchi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if printbool==1
    fprintf('\n------------------------------------------ Linear Fitting ------------------------------------------\nDegrees of Freedom=%i',(length(x)-N))
end


N2   = length(x);                  % Groesse des Datensatzes.
n = length(funccell);
doF = N2-n-npolyfit;                        % Degrees of Freedom.

%%% INITALISIEREN VON cerr, G und G2 %%%

cerr = zeros(n,1);
G    = zeros(n,N2);
G2   = zeros(n,N2);

%%% BERECHNEN VON g und g2 %%%

for j=0:npolyfit-1
    funcstr2 = sprintf('funccell{%i+1+n}= @(x) x.^(%i);' ,[j j]);
    eval(funcstr2);
end

n2 = n+npolyfit;

for alpha=1:n2,
    func = funccell{alpha};
    G(alpha,:)  = func(x);    % Siehe Gl. (2).
    G2(alpha,:) = G(alpha,:)./(yerr.^2); % Siehe Gl. (3).
end

%%%%%%%%% Test whether functions are linear independent %%%%%%%%%
try
    if rank(G)<n2
        Gtest = G(1:npolyfit+n,1:npolyfit+n)';
        [rub1 geig V] = svd(Gtest);
        rankG=rank(G);
        for i=rankG+1:(npolyfit+n)
            lincofm(:,i-rankG)=V(:,i);
        end
        ndep = (npolyfit+n)-rankG;
        fprintf('############################## Warning: ########################### \n There are %i linear dependencies in your functions. This leads to bad results.\n The linear dependencies are:\n\n ',ndep)
        for i=1:ndep
            lincof=lincofm(:,i);
            d1 = abs(lincof)>0.00001*max(abs(lincof));
            d2 = [1:length(d1)]'.*d1;
            linfunc = d2(d2~=0);
            lincof2=lincof(d1);
            lincof2=lincof2/lincof2(1);
            for j=1:length(linfunc)
                if lincof2(j)>0
                    fprintf('+')
                end
                fprintf('%2.3f*',lincof2(j))
                funcstr2 = func2str(funccell{linfunc(j)});
                fprintf('%s ',funcstr2(5:length(funcstr2)))
                if j==length(linfunc)
                    fprintf('= 0\n')
                end
            end

            if i~=(npolyfit+n)-rankG
                fprintf('\n')
            end
        end
        fprintf('\n#############################################################\n')
    end
catch ME
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BERECHNEN VON A und b %%%

A = G*G2';                     % Siehe Gl. (4).
b = G2*y';                     % Siehe Gl. (5).

Ainv = A\eye(n2);               % Einmal invertieren...

%%% BERECHNEN DER FITPARAMETER %%%

c = Ainv*b;                       % Fitparameter, siehe Gl. (6).

Kov = Ainv;                       % Kovarianz-Matrix, Gl. (7), Kov <-> M.

for k=1:n+npolyfit                         % Siehe Gl. (9.20).
    cerr(k,1) = sqrt(Kov(k,k));
end

%%%%%%%%%%%% Uebergabe an Franziska Flegels Programm %%%%%%%%%%%%

for i=1:n2
    eval(['func' num2str(i) '=funccell{' num2str(i) '};'])
end

funcoutstr = 'funcout = @(x)';
for i=1:n2
    funcoutstr = [funcoutstr 'c(' num2str(i) ')*func' num2str(i) '(x)' ] ;
    if i<n2
        funcoutstr = [funcoutstr '+'];
    else
        funcoutstr = [funcoutstr ';'];
    end
end


eval(funcoutstr);
chimin = chi2(x,y,yerr,funcout);
beta = c';
betaerr = cerr';
chi2dof=chimin/doF;
Q=1 -gammainc(chimin/2,doF/2);

if printbool
    fprintf('    chi^2/dof = %1.2f\n',chi2dof)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





if any(isnan(beta)) || any(isnan(betaerr)) || any(isinf(beta))
    error('No results could be found. This is probably due to a bad function input.\n Test whether the values of the input functions are not near zero')
end

% Plotting


if plotbool==1;

    dtextint = get(0, 'defaultTextInterpreter');
    set(0, 'defaultTextInterpreter', 'latex');

    if ~isempty(font)
        dfonta=get(0,'defaultAxesFontName');
        dfontt=get(0,'defaultTextFontName');
        set(0,'defaultAxesFontName', font)
        set(0,'defaultTextFontName', font)
    end

    if isempty(stringcell)
        stringcell={'x-axis' 'y-axis'};
        for i=1:N
            stringcell{i+2}=sprintf('$c_{%i}$',i);
        end
    end

    N2=length(stringcell)    ;
    if 2<=N2<N+2
        for i=1:N+2-N2
            stringcell{i+N2}=sprintf('$c_{%i}$',i);
        end
    end


end
% Residuen
residue=(y-funcout(x))./yerr;
mresidue=mean(residue);


% Normalverteilungstests
if length(residue)>10 && printbool && exist('kstest','file') && normalitytestbool
    [kstestbin pks]=kstest(residue,[],0.05,'unequal');
    titlestr='\nResidues:\n';
    if kstestbin
        titlestr=[titlestr 'Residues are not standard normally distributed with 5%% Significance. p=' num2str(pks,3) '\n' ];
    else
        titlestr=[titlestr 'The null hypothesis that the sample came from a standard normally distributed \npopulation could not be rejected with 5%% significance. p=' num2str(pks,3) '\n' ];
    end

    if exist('swtest','file')

        [normtest p]=swtest(residue,0.05,0);
        titlestr=[titlestr '\n'];
        if normtest
            titlestr=[titlestr 'Residues are not normally distributed with 5%% Significance. p=' num2str(p,3) '\n' ];
        else
            titlestr=[titlestr 'The null hypothesis that the sample came from a normally distributed \npopulation could not be rejected with 5%% significance. p=' num2str(p,3) '\n' ];
        end


    end

    fprintf(titlestr)
end
if plotbool

    x=x*rescalex;
    y=y*rescaley;
    yerr=yerr*rescaley;

    breiterr=(x(length(x))-x(1))/100;
    figure(2);
    clf;
    set(gca,'FontSize',12)
    plot(x,mresidue*ones(size(x)),'r','linewidth',2);
    hold on
    resfig=plot(x,residue,'ko');
    hold off


    axisx1=(min(x) -(max(x)-min(x))/30);
    axisx2=(max(x) +(max(x)-min(x))/30);

    axisy1=min(residue) -(max(residue)-min(residue))/40;
    axisy2=max(residue) +(max(residue)-min(residue))/40;

    axis([axisx1 axisx2 axisy1 axisy2])


    xlabel(stringcell{1},'fontsize',16)
    h2=ylabel('Residues','fontsize',16);
    set(h2,'unit','character')
    legend('Mean of residues','Residues','Location','NorthOutside')
    set(resfig                            , ...
        'LineWidth'       , 1           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 5           , ...
        'MarkerEdgeColor' , [.3 .3 .3]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , ygridstr  , ...
        'XGrid'       , xgridstr  , ...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'LineWidth'   , 1         );

    set(gcf, 'PaperPositionMode', 'manual');    % Macht, dass du das selbst einstellen darfst.
    set(gcf, 'PaperUnits', 'centimeters');                 % Macht, dass du die Groessen in Inches angeben kannst. (Es geht auch z.B. 'centimeters')
    set(gcf, 'PaperSize', format);                       % Einstellen der gewuenschten Groesse.
    set(gcf, 'PaperPosition', [0 0 format(1)*1.05 format(2)]);

    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',[-8 13])


    % Fit mit Daten
    figure(1);
    clf;
    set(gca,'FontSize',12)

    if isempty(axisin)

        axisx1=(min(x) -(max(x)-min(x))/30);
        axisx2=(max(x) +(max(x)-min(x))/30);

        [miny Iminy]=min(y);
        [maxy Imaxy]=max(y);

        axisy1=min(y)-yerr(Iminy) -(max(y)-min(y))/40;
        axisy2=max(y)+yerr(Imaxy) +(max(y)-min(y))/40;

        axisin=[axisx1 axisx2 axisy1 axisy2];
    elseif (isnumeric(axisin) && length(axisin)==4)
        axisx1=axisin(1);
        axisx2=axisin(2);
    else
        warning('Matlab:wnonlinfit:badaxisin','Bad axis input. Default is used')

        axisx1=(min(x) -(max(x)-min(x))/30);
        axisx2=(max(x) +(max(x)-min(x))/30);

        [miny Iminy]=min(y);
        [maxy Imaxy]=max(y);

        axisy1=min(y)-yerr(Iminy) -(max(y)-min(y))/40;
        axisy2=max(y)+yerr(Imaxy) +(max(y)-min(y))/40;

        axisin=[axisx1 axisx2 axisy1 axisy2];
    end

    xplot=linspace(axisx1,axisx2,500)';
    endfitplot=funcout(xplot/rescalex)*rescaley;

    if isempty(plotmodestr)
        hFit=plot(xplot,endfitplot,'k');
    elseif strcmpi(plotmodestr,'semilogx')
        hFit=semilogx(xplot,endfitplot,'k');
    elseif strcmpi(plotmodestr,'semilogy')
        hFit=semilogy(xplot,endfitplot,'k');
    elseif strcmpi(plotmodestr,'loglog')
        hFit=loglog(xplot,endfitplot,'k');
    else
        fprintf('invalid input for plotmode. Default used.')
        hFit=plot(xplot,endfitplot,'k');
    end
    axis(axisin);

    % Text positioning
    max1=max(endfitplot(1:250));
    max2=max(endfitplot(251:498));
    if isempty(position)  && max1<=max2
        position(1)=axisin(1)+(x(length(x))-x(1))/40;
        position(2)=axisin(4)-(axisin(4)-axisin(3))/40;
    elseif isempty(position) && max1>max2
        position(1)=axisin(2)-(axisin(2)-axisin(1))/2;
        position(2)=axisin(4)-(axisin(4)-axisin(3))/40;
    end


    hold on
    hE=errorbar(x,y,yerr ,'bo','MarkerSize',5);
    hold off
    xlabel(stringcell{1},'Interpreter','Latex','fontsize',16);
    h1=ylabel(stringcell{2},'Interpreter','Latex','fontsize',16);

    % Positioning of x/y labels    
    
    set(h1,'unit','character')
    set(h1,'Position',get(h1,'Position') +  [-3 0 0])
    
    
    strings=cell(1,N);
    for i=1:N
        strings{i}=stringcell{i+2};
    end

try
    parameters=paramtext(textheader,beta,betaerr,errprec,chi2dof,Q,strings);
catch
    error('Function paramtext is missing. Please add into matlab or local folder')
end

    if ~chibool
        parameters=parameters(1:length(parameters)-3);
    end


   htext =  text(position(1),position(2),parameters,'VerticalAlignment',...
        'top','fontsize',14,'BackgroundColor',[0.95 0.95 0.95],'EdgeColor','k','Margin',1);



    set(hFit,'LineWidth',2);
    if isempty(plotmodestr)
        hE_c                   = ...
            get(hE     , 'Children'    );
        errorbarXData          = ...
            get(hE_c(2), 'XData'       );
        errorbarXData(4:9:end) = ...
            errorbarXData(1:9:end) - breiterr;
        errorbarXData(7:9:end) = ....
            errorbarXData(1:9:end) - breiterr;
        errorbarXData(5:9:end) = ...
            errorbarXData(1:9:end) + breiterr;
        errorbarXData(8:9:end) = ...
            errorbarXData(1:9:end) + breiterr;
        set(hE_c(2), 'XData', errorbarXData);
    end
    set(hE                            , ...
        'LineStyle'       , 'none'      , ...
        'Marker'          , '.'         , ...
        'Color'           , [0.3 0.3 0.3]  );

    set(hE                            , ...
        'LineWidth'       , 1           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 5           , ...
        'MarkerEdgeColor' , [.3 .3 .3]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    set( gca                       , ...
        'FontName'   , 'Helvetica' );

    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , ygridstr  , ...
        'XGrid'       , xgridstr  , ...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'LineWidth'   , 1         );


    set(gcf, 'PaperPositionMode', 'manual');    % Macht, dass du das selbst einstellen darfst.
    set(gcf, 'PaperUnits', 'centimeters');                 % Macht, dass du die Groessen in cm angeben kannst. (Es geht auch z.B. 'centimeters')
    set(gcf, 'PaperSize', format);                       % Einstellen der gewuenschten Groesse.
    set(gcf, 'PaperPosition', [0 0 format(1)*1.05 format(2)]);


    if ~isempty(contourcell)
        contourfunc(Kov,contourcell,beta,betaerr,format,strings)
    end

    set(0, 'defaultTextInterpreter', dtextint);
    if ~isempty(font)
        set(0,'defaultAxesFontName', dfonta)
        set(0,'defaultTextFontName', dfontt)
    end
figure(1);
end
beta=beta';
betaerr=betaerr';



end

%%%%%% Other functions

function [x] = logicalvar(x,varargin )
ww=varargin{1};

if isnumeric(x)
    if x==1
        x=true;
    elseif x==0
        x=false;
    else
        warning('MATLAB:locicalvar:inaprinput',['Inapropriate input for ',ww,'. Standard value true is used'])
        x=true;
    end
end

if ischar(x)
    if any(strcmpi(x,{'y','yes','ja','oui','si','on'}))
        x=true;
    elseif any(strcmpi(x,{'n','no','nein','non','off'}))
        x=false;
    else
        warning('MATLAB:locicalvar:inaprinput',['Inapropriate input for ',ww,'. Standard value true is used'])
        x=true;
    end
end

end


function [chi] = chi2(x,y,err,func)

% Returns chi^2 value for chi=chi2(x,y,yerr,function)

if length(err)==1
    err=err*ones(size(y));
end
if length(x)~=length(y) || length(y)~=length(err)
    error('Data must have same size')
end
chi=sum( ((func(x)-y)./err).^2) ;
end

function contourfunc(Kovm,contourcell,betam,betaerrm,format,strings)
n = length(contourcell);
for i=1:n

    var = contourcell{i};
    Kov = [Kovm(var(1),var(1)) Kovm(var(1),var(2));Kovm(var(1),var(2)) Kovm(var(2),var(2))] ;
    beta = [betam(var(1));betam(var(2))];
    betaerr = [betaerrm(var(1));betaerrm(var(2))];

    [A, B] = eig(Kov^-1);
    t = 0:0.01:2*pi;
    figure(i+2)

    sigma1 = repmat(beta,1,length(t))+A * [1/sqrt(B(1,1))*cos(t);1/sqrt(B(2,2))*sin(t)];
    sigma2 = repmat(beta,1,length(t))+A * 2*[1/sqrt(B(1,1))*cos(t);1/sqrt(B(2,2))*sin(t)];
    sigma3 = repmat(beta,1,length(t))+A * 3*[1/sqrt(B(1,1))*cos(t);1/sqrt(B(2,2))*sin(t)];


    plot(beta(1),beta(2),'k.','linewidth',2);
    hold on
    plot(sigma1(1,:),sigma1(2,:),'linewidth',2,'color',[0 0 0]);
    plot(sigma2(1,:),sigma2(2,:),'linewidth',2,'color',[0.5 0.5 0.5]);
    plot(sigma3(1,:),sigma3(2,:),'linewidth',2,'color',[0.8 0.8 0.8]);
    hold off
    axis([beta(1)-3.3*betaerr(1) beta(1)+3.3*betaerr(1) beta(2)-3.3*betaerr(2) beta(2)+3.3*betaerr(2)])

    
    h_legend = legend(' Expectation value',' 1 $\sigma$',' 2 $\sigma$', ' 3 $\sigma$');

    legend(h_legend,'boxoff')
    set(h_legend,'FontSize',14);
    set(h_legend, 'Interpreter', 'latex')
    if Kov(1,2)<0
        set(h_legend,'Location', 'NorthEast')
    else
        set(h_legend, 'Location', 'NorthWest')
    end

    xlabel(['Fitparameter ' strings{var(1)}],'Interpreter','Latex','fontsize',16);
    ylabel(['Fitparameter ' strings{var(2)}],'Interpreter','Latex','fontsize',16);

    set(gca, ...
        'FontSize'    ,  12      , ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'LineWidth'   , 1         );


    title(['Contour plot: Bivariate pdf of ' strings{var(1)} ' and ' strings{var(2)}],'Interpreter','Latex','fontsize',16)


    set(gcf, 'PaperPositionMode', 'manual');    % Macht, dass du das selbst einstellen darfst.
    set(gcf, 'PaperUnits', 'centimeters');                 % Macht, dass du die Groessen in cm angeben kannst. (Es geht auch z.B. 'centimeters')
    set(gcf, 'PaperSize', format);                       % Einstellen der gewuenschten Groesse.
    set(gcf, 'PaperPosition', [0 0 format(1)*1.05 format(2)]);

end
end