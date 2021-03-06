% Weighted nonlinear fit by chi^2 minimization.
%
% Last update: 12.07.13
%
% wnonlinfit returns parameters,errors,chi2dof,probablility and a vector
% containing all chi^2 that accured during the search. If the chi^2/dof
% values are worse than a defined value (default is 2) random changes are
% done to the initial guesses in order to find a better guess.
%
%   [parameters errors chi2dof]=wnonlinfit(x,y,yerr,func,beta0,options)
%
%
% x,y,yerr are the data variables.
%
% func is the model that is to be used. It must be a function handle in the
% form func=@(xdata,betav) ... . So that the result of func(x,beta0) is a vector of size(x).
% betav are the variables that are to be fitted.
%
% beta0 is the initial guess.
%
% Important: All text input is in latex
%
% Options:
%
% Options are entered like: wnonlinfit(x,y,yerr,func,beta0,'optionname',value,'optionname2',value2)
%
%
%
% chitol: scalar in ]0,inf[                     Defines chi^2/dof tolerance.
%                                               Routine searches for better guesses until chi^2<chitol*dof
%                                               Default: 2
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
% 
% normalitytests:                               Defines whether
% 'on' , true , 1 or 'off' , false , 0          normalitytests are printed
%                                               on the screen
% 
% 
% 
% contourplt: matrices of 1x2 numbers			For every matrix a correlation
%												contourplot between the defined variables is created 
%												Default: no contourplot
%												Example: {[1 2],[1 3]}
%												creates correlation plots
%												of variable (1 and 2) and (1
%												and 3)
% 
% Written by: Jannick Weisshaupt (jannick@physik.hu-berlin.de)
% Paramtext by: Franziska Flegel
%
% Feel free to share, change or do whatever you want with this script.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%% Example (should run on its own. Just copy to test and understand):
%
% clear;
%
% % you need to specify the function you want to use for fitting
% % the syntax is
% % fitfunc(vector of xvalues,vector of fitparameters) = @(x,beta) ...
% fitfunc=@(x,betav) betav(1)*x.^3.*1./(exp(x.*betav(2))-1);
%
% % you also need a first guess for the fitparameters:
% beta0=[3 2];
%
% %%%% Creation of data %%%%%%
% beta1=[1.5 0.75];
% error=0.1 ;
% x=[10*eps:0.5:20];
% y=fitfunc(x,beta1);
% y=y+(randn(size(y)))*error;
% yerr=error.*ones(size(y));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % normally the input data has to be vector of xvalues, vector of yvalues,
% % vector of errors. If you don't have errors you can just put a 1 (or
% % arbitrary value) instead of a vector
%
%
%
% %%%%% Optional input %%%%%
% textposition=[];    % define position of legend
%
% % define header of legend
% headercell={'Nonlinear Fit'  '$f(x)=A_0 \cdot \frac{1}{\exp{(x \beta)}-1}$' ''};
% % define axis labels and fitparameter names
% mylabel={'xaxis [eV]','yaxis [Counts]','$A_0$','$\beta$'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% % actual execution of wnonlinfit. See "help wnonlinfit" for more options.
% [beta betaerr chi prob chiminvec]=wnonlinfit(x,y,yerr,fitfunc,beta0...          % main arguments
%     ,'label',mylabel,'position',textposition,'header',headercell,'errprec',2);  % optional arguments
%
%
% print -dpdf -cmyk test.pdf
%
% figure(2)
% print -dpdf -cmyk testresidue.pdf





function [beta,betaerr,Kov,chi2dof,prob,chiminvec] = wnonlinfit(x,y,yerr,func,beta0,varargin)
N=length(beta0);

%%%%%%%%%%%%%%%%%%%%%%% Testing %%%%%%%%%%%%%%%%%%%
ytest=func(x,beta0);

if any(~isnumeric(x)) || any(~isnumeric(y)) || any(~isnumeric(yerr))
    error('There are non numerical entries in your data (x,y or yerr)')
end

if length(yerr)==1
    yerr=yerr *ones(size(y));
end

if length(x)~=length(y) || length(y)~=length(yerr)
    error('Data must have same size')
end

if any(size(y)~=size(ytest))
    error('Function values must have same size as y values')
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
    warning('MATLAB:wnonlinfit:small_error_warning','Fool! Errors are smaller than machine accuracy.\nThis might lead to very strange results')
end

if any(isnan(ytest)) || any(isinf(ytest))
    error('MATLAB:wnonlinfit:badguess','Test your function! Infinite or NAN values encountered with your starting guesses.')
end


% Optional arguments

p = inputParser;

p.addParamValue('plot',true);
p.addParamValue('label',[]);
p.addParamValue('chitol',2)
p.addParamValue('position',[])
p.addParamValue('print',1)
p.addParamValue('format',[16 14])
p.addParamValue('grid',false)
p.addParamValue('header','Nonlinear Fit')
p.addParamValue('printchi',true)
p.addParamValue('axis',[]);
p.addParamValue('errprec',2);
p.addParamValue('rescalex',1,@(x) isnumeric(x) && length(x)==1);
p.addParamValue('rescaley',1,@(x) isnumeric(x) && length(x)==1);
p.addParamValue('plotmode',[]);
p.addParamValue('font',[]);
p.addParamValue('contourplot',{});
p.addParamValue('normalitytests',false);


p.parse(varargin{:});

stringcell=p.Results.label;
position=p.Results.position;
chitol=p.Results.chitol;
plotbool=p.Results.plot;
printbool=p.Results.print;
gridbool=p.Results.grid;
format=p.Results.format;
textheader=p.Results.header;
chibool=p.Results.printchi;
axisin=p.Results.axis;
errprec=p.Results.errprec;
font=p.Results.font;
rescalex=p.Results.rescalex;
rescaley=p.Results.rescaley;
plotmodestr=p.Results.plotmode;
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

plotbool=logicalvar(plotbool,'plot');
printbool=logicalvar(printbool,'print');
chibool=logicalvar(chibool,'printchi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printbool==1
    fprintf('\n------------------------------------------ Nonlinear Fitting ------------------------------------------\nDegrees of Freedom=%i',(length(x)-N))
end
i=1;
testbin=0;
dof=(length(y)-N);

if printbool==1
    fprintf('     Iteration:   ')
    fprintf('%i',i)
end
options = optimset('Display','off');


if exist('lsqnonlin','file')
    [betavec(1,:), chiminvec(1,1)]=lsqnonlin(@(betav) (y-func(x,betav))./yerr ,beta0,[],[],options);
else
    [betavec(i,:), chiminvec(i,1)]=fminsearch(@(betav) chi2(x,y,yerr,@(y) func(y,betav)),beta0.*(1+((i-51)/100+1)*randn(size(beta0))),options);

end

if chiminvec(1)>dof*chitol
    for i=2:50
        if printbool==1;
            if i<10
                fprintf('\b');
            else
                fprintf('\b\b');
            end
            fprintf('%i',i)
        end

        if exist('lsqnonlin','file')
            [betavec(i,:), chiminvec(i,1)]=lsqnonlin(@(betav) (y-func(x,betav))./yerr ,beta0.*(1+randn(size(beta0))/5),[],[],options);
        else
            [betavec(i,:), chiminvec(i,1)]=fminsearch(@(betav) chi2(x,y,yerr,@(y) func(y,betav)),beta0.*(1+((i-51)/100+1)*randn(size(beta0))),options);
        end



        if chiminvec(i)<dof*chitol+sum(y)*eps
            if printbool==1
                fprintf('     Success (Chi2dof<%0.1f) ',chitol)
            end
            break
        end
        if i==10 && max(func(x,betavec(i,:)))-min(func(x,betavec(i,:)))>10*eps*min(func(x,betavec(i,:)))
            if numel(find(abs((min(chiminvec)-chiminvec)/chiminvec(1))<0.01))>=(numel(chiminvec))

                if printbool==1
                    fprintf('     Success (Same chi^2 (0.1%% Tol.) for 10 iterations)')
                end

                testbin=1;
                break
            end
        end
        if i==22 && max(func(x,betavec(i,:)))-min(func(x,betavec(i,:)))>10*eps*min(func(x,betavec(i,:)))
            if numel(find(abs((min(chiminvec)-chiminvec)/chiminvec(1))<0.05))>=(numel(chiminvec)-2)
                if printbool==1
                    fprintf('     Success (Same chi^2 (1%% Tol.) for 20 iterations)')
                end
                testbin=1;
                break
            end
        end
    end
elseif printbool==1
    fprintf('     Success (Chi2dof<%0.1f)',chitol)
end

[chimin, nchimin]=min(chiminvec);
beta=betavec(nchimin,:);

if chimin>dof*chitol+sum(y)*eps && testbin==0
    for i=51:250
        if printbool==1
            if i<100
                fprintf('\b\b');
            else
                fprintf('\b\b\b');
            end
        end
        if printbool==1
            fprintf('%i',i)
        end
        if exist('lsqnonlin','file')
            [betavec(i,:), chiminvec(i,1)]=lsqnonlin(@(betav) (y-func(x,betav))./yerr ,beta0.*(1+((i-51)/100+1)*randn(size(beta0))),[],[],options);
        else
            [betavec(i,:), chiminvec(i,1)]=fminsearch(@(betav) chi2(x,y,yerr,@(y) func(y,betav)),beta0.*(1+((i-51)/100+1)*randn(size(beta0))),options);
        end
        if chiminvec(i,1)<((i-51)/10+1)*dof*chitol
            if printbool==1
                fprintf('     Success (Chi2dof<%0.2f)',((i-51)/10+1)*chitol)
            end
            break
        end

        if i==100 && max(func(x,betavec(i,:)))-min(func(x,betavec(i,:)))>10*eps*min(func(x,betavec(i,:)))
            if numel(find(abs((min(chiminvec)-chiminvec)/min(chiminvec))<0.05))>=(numel(chiminvec)-10)
                fprintf('     Success (Same chi^2 (5%% Tol) for 90 iterations)')
                break
            end
        end

        if i==250
            fprintf('    Did not fulfill criterias before reaching the maximum iteration number')
        end

    end

    [chimin, nchimin]=min(chiminvec);
    beta=betavec(nchimin,:);
end


chiminvec=chiminvec/dof;
chi2dof=chimin/dof;
Q = 1 - gammainc(chimin/2,dof/2);

if printbool==1
    fprintf('   Chi2dof=%1.2f  Q=%1.2f\n',[chi2dof Q])
end

[betaerr, Kov]=errorfit(x,y,yerr,func,beta);

if any(isnan(betaerr)) && printbool==1
    warning('MATLAB:wnonlinfit:fzeronan','Fzero failed at finding errors. Test your function!')
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
            stringcell{i+2}=sprintf('$c_%i$',i);
        end
    end

    N2=length(stringcell)    ;
    if 2<=N2<N+2
        for i=1:N+2-N2
            stringcell{i+N2}=sprintf('$c_%i$',i);
        end
    end


end
% Residuen
residue=(y-func(x,beta))./yerr;
mresidue=mean(residue);


% Normalverteilungstests
if length(residue)>10 && printbool && exist('kstest','file') && normalitytestbool
    [kstestbin, pks]=kstest(residue,[],0.05,'unequal');
    titlestr='\nResidues:\n';
    if kstestbin
        titlestr=[titlestr 'Residues are not standard normally distributed with 5%% Significance. p=' num2str(pks,3) '\n' ];
    else
        titlestr=[titlestr 'The null hypothesis that the sample came from a standard normally distributed \npopulation could not be rejected with 5%% significance. p=' num2str(pks,3) '\n' ];
    end

    if exist('swtest','file')

        [normtest, p]=swtest(residue,0.05,0);
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

        [miny, Iminy]=min(y);
        [maxy, Imaxy]=max(y);

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
    endfitplot=func(xplot/rescalex,beta)*rescaley;

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
    h1 = ylabel(stringcell{2},'Interpreter','Latex','fontsize',16);
    % Positioning of x/y labels    
    
    set(h1,'unit','character')
    set(h1,'Position',get(h1,'Position') +  [-3 0 0])
   
   
    strings=cell(1,N);
    for i=1:N
        strings{i}=stringcell{i+2};
    end


    parameters=paramtext(textheader,beta,betaerr,errprec,chi2dof,Q,strings);
    if ~chibool
        parameters=parameters(1:length(parameters)-3);
    end


    text(position(1),position(2),parameters,'VerticalAlignment',...
        'top','fontsize',14,'BackgroundColor',[0.95 0.95 0.95],'EdgeColor','k','Margin',1)



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


    set(0, 'defaultTextInterpreter', dtextint);
    if ~isempty(font)
        set(0,'defaultAxesFontName', dfonta)
        set(0,'defaultTextFontName', dfontt)
    end

    if ~isempty(contourcell)
        contourfunc(Kov,contourcell,beta,betaerr,format,strings)
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


function [betaerr, Kov] = errorfit(x,y,yerr,func,beta,chimin)
if nargin>5

    chimin2=chi2(x,y,yerr,@(y) func(y,beta));

    if abs((chimin2-chimin)/chimin)>0.001
        warning('MATLAB:errorfit:wrongchimin','wrong chimin');
    end

end

if nargin<6
    chimin=chi2(x,y,yerr,@(y) func(y,beta));
end

[a1 b1]=size(x);
if a1<b1,
   x=x';
   y=y';
   yerr=yerr';
end

N=length(beta);
n = length(x);
dof = n-N;

funcvec = @(beta) func(x,beta);
F = jacobianest(funcvec,beta);
G=zeros(n,N);
for i=1:N
    G(:,i) = F(:,i)./(yerr.^2);
end

Kov = (F'*G)^-1;
% ssquare = dof^-1 *yerr'*(eye(n,n)-F*C*F')*yerr;
% Kov = 1.405^2*ssquare*C;
betaerr = sqrt(diag(Kov))';

end

function [jac,err] = jacobianest(fun,x0)
% gradest: estimate of the Jacobian matrix of a vector valued function of n variables
% usage: [jac,err] = jacobianest(fun,x0)
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/6/2007

% get the length of x0 for the size of jac
nx = numel(x0);

MaxStep = 100;
StepRatio = 2.0000001;

% was a string supplied?
if ischar(fun)
    fun = str2func(fun);
end

% get fun at the center point
f0 = fun(x0);
f0 = f0(:);
n = length(f0);
if n==0
    % empty begets empty
    jac = zeros(0,nx);
    err = jac;
    return
end

relativedelta = MaxStep*StepRatio .^(0:-1:-25);
nsteps = length(relativedelta);

% total number of derivatives we will need to take
jac = zeros(n,nx);
err = jac;
for i = 1:nx
    x0_i = x0(i);
    if x0_i ~= 0
        delta = x0_i*relativedelta;
    else
        delta = relativedelta;
    end

    % evaluate at each step, centered around x0_i
    % difference to give a second order estimate
    fdel = zeros(n,nsteps);
    for j = 1:nsteps
        fdif = fun(swapelement(x0,i,x0_i + delta(j))) - ...
            fun(swapelement(x0,i,x0_i - delta(j)));

        fdel(:,j) = fdif(:);
    end

    % these are pure second order estimates of the
    % first derivative, for each trial delta.
    derest = fdel.*repmat(0.5 ./ delta,n,1);

    % The error term on these estimates has a second order
    % component, but also some 4th and 6th order terms in it.
    % Use Romberg exrapolation to improve the estimates to
    % 6th order, as well as to provide the error estimate.

    % loop here, as rombextrap coupled with the trimming
    % will get complicated otherwise.
    for j = 1:n
        [der_romb,errest] = rombextrap(StepRatio,derest(j,:),[2 4]);

        % trim off 3 estimates at each end of the scale
        nest = length(der_romb);
        trim = [1:3, nest+(-2:0)];
        [der_romb,tags] = sort(der_romb);
        der_romb(trim) = [];
        tags(trim) = [];

        errest = errest(tags);

        % now pick the estimate with the lowest predicted error
        [err(j,i),ind] = min(errest);
        jac(j,i) = der_romb(ind);
    end
end

end % mainline function end

% =======================================
%      sub-functions
% =======================================
function vec = swapelement(vec,ind,val)
% swaps val as element ind, into the vector vec
vec(ind) = val;

end % sub-function end

% ============================================
% subfunction - romberg extrapolation
% ============================================
function [der_romb,errest] = rombextrap(StepRatio,der_init,rombexpon)
% do romberg extrapolation for each estimate
%
%  StepRatio - Ratio decrease in step
%  der_init - initial derivative estimates
%  rombexpon - higher order terms to cancel using the romberg step
%
%  der_romb - derivative estimates returned
%  errest - error estimates
%  amp - noise amplification factor due to the romberg step

srinv = 1/StepRatio;

% do nothing if no romberg terms
nexpon = length(rombexpon);
rmat = ones(nexpon+2,nexpon+1);
% two romberg terms
rmat(2,2:3) = srinv.^rombexpon;
rmat(3,2:3) = srinv.^(2*rombexpon);
rmat(4,2:3) = srinv.^(3*rombexpon);

% qr factorization used for the extrapolation as well
% as the uncertainty estimates
[qromb,rromb] = qr(rmat,0);

% the noise amplification is further amplified by the Romberg step.
% amp = cond(rromb);

% this does the extrapolation to a zero step size.
ne = length(der_init);
rhs = vec2mat(der_init,nexpon+2,ne - (nexpon+2));
rombcoefs = rromb\(qromb'*rhs);
der_romb = rombcoefs(1,:)';

% uncertainty estimate of derivative prediction
s = sqrt(sum((rhs - rmat*rombcoefs).^2,1));
rinv = rromb\eye(nexpon+1);
cov1 = sum(rinv.^2,2); % 1 spare dof
errest = s'*12.7062047361747*sqrt(cov1(1));

end % rombextrap


% ============================================
% subfunction - vec2mat
% ============================================
function mat = vec2mat(vec,n,m)
% forms the matrix M, such that M(i,j) = vec(i+j-1)
[i,j] = ndgrid(1:n,0:m-1);
ind = i+j;
mat = vec(ind);
if n==1
    mat = mat';
end

end % vec2mat

function contourfunc(Kovm,contourcell,betam,betaerrm,format,strings)
n = length(contourcell);
for i=1:n

    var = contourcell{i};
    Kov = [Kovm(var(1),var(1)) Kovm(var(1),var(2));Kovm(var(1),var(2)) Kovm(var(2),var(2))] ;
    beta = [betam(var(1));betam(var(2))];
    betaerr = [betaerrm(var(1));betaerrm(var(2))];

    [A, B] = eig(Kov^-1);
    t = 0:0.001:2*pi;
    figure(i+2);

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

    %         ylabh = get(gca,'YLabel');
    %     set(ylabh,'Position',[13 13])


    title(['Contour plot: Bivariate pdf of ' strings{var(1)} ' and ' strings{var(2)}],'Interpreter','Latex','fontsize',16)


    set(gcf, 'PaperPositionMode', 'manual');    % Macht, dass du das selbst einstellen darfst.
    set(gcf, 'PaperUnits', 'centimeters');                 % Macht, dass du die Groessen in cm angeben kannst. (Es geht auch z.B. 'centimeters')
    set(gcf, 'PaperSize', format);                       % Einstellen der gewuenschten Groesse.
    set(gcf, 'PaperPosition', [0 0 format(1)*1.05 format(2)]);

end


end