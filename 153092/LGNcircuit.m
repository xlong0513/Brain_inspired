function LGNcircuit

% LGN_GUI Give input vector/function to the LGN relay cell and choose
% parameters. The GUI plots the resulting LGN response including separate
% components of the response, if requested.
    
    % including the data files and data file path
    file_reg = ['dataFiles' filesep '*.dat'];
    file_struct = dir(file_reg);
    filenames = cell(1, length(file_struct));
    for ifile = 1:length(file_struct)
        filenames(1,ifile)=cellstr(file_struct(ifile).name);
    end
    
    data_path = ['dataFiles' filesep];
    addpath(data_path);
    
    % Setting default parameter values    
    lgn_struct=struct('eta_ffi',0.7,'tau_rg',20,'tau_rig',50,'Delta_rig',0,...
        'l_rON',-1,'w_fbON',0,'w_fbOFFx',0,'tau_rc',20,'Delta_fb',10,...
        'l_cON',0.1,'l_cOFF',-0.1);
    in_struct=struct('form','vec',...
        'interpol_method','linear',...
        'data_file',' ','mode','gui');
    
    % Set properties for the parameters table
    rowname = {'eta_ffi', 'tau_rg', 'tau_rig', 'Delta_rig', 'l_rON' 'w_fbON',...
        'w_fbOFFx', 'tau_rc', 'Delta_fb', 'l_cON', 'l_cOFF'};
    rowname_tab = {'eta_ffi', 'tau_rg (ms)', 'tau_rig (ms)', 'Delta_rig (ms)',...
        'l_rON','w_fbON','w_fbOFFx', 'tau_rc (ms)', 'Delta_fb (ms)', 'l_cON',...
        'l_cOFF'};
    
    dat = {lgn_struct.eta_ffi;...
        lgn_struct.tau_rg;...
        lgn_struct.tau_rig;...
        lgn_struct.Delta_rig;...
        lgn_struct.l_rON;...
        lgn_struct.w_fbON;...
        lgn_struct.w_fbOFFx;...
        lgn_struct.tau_rc;...
        lgn_struct.Delta_fb;...
        lgn_struct.l_cON;...
        lgn_struct.l_cOFF};
    columnname = {'Value'};
    columnformat = {'numeric'};
    columneditable = true;
    columnwidth = {80};
    
    ih1 = 1;
    
    % Initialize and hide GUI at construction.
    f_width=900;
    f_height=900;
    f = figure('Visible','off','Tag','gui',...
        'Position',[0,100,f_width,f_height]);
    %    'BackgroundColor',[0.5 0.5 0.5],...
    % Define an axes that covers the whole figure
    %ax_root = axes('position', [0,0,1.0,1.0]); 
    
    % Construct the components
    bcgp = uipanel('Parent',f,'Units','pixels',...
        'BackgroundColor',[.7 .7 .7],...
        'Position',[0,0,f_width,f_height]); 
    mp = uipanel('Parent',f,'Title','LGNcircuit','Units','pixels',...
        'Position',[570,160,315,720]);
%     inputform_text = uicontrol(mp,'Units','pixels','Style','text',...
%         'String','Select input format','Position',[10,730,160,20]);
%     inputform_pop = uicontrol(mp,'Units','pixels','Style','popupmenu',...
%         'String',{'Data file'},...
%         'Position',[10,695,180,35],...
%         'BackgroundColor','w',...
%         'Callback',{@popup_callback});
%     ipf_help = uicontrol(mp,'Units','pixels','Style','pushbutton',...
%         'String','Help',...
%         'Position',[200,704,105,30],...
%         'Tag','ipf',...
%         'Callback',{@helpbutton_callback});
    inputdata_text = uicontrol(mp,'Units','pixels','Style','text',...
        'String','Data file',...
        'Position',[10,665,100,20]);
    inputdata_pop = uicontrol(mp,'Units','pixels',...
        'Style','popupmenu',...
        'String',filenames,...
        'BackgroundColor','w',...
        'Position',[10,630,180,30],...
        'Callback',{@datapop_callback});
    ipd_help = uicontrol(mp,'Units','pixels','Style','pushbutton',...
        'String','Help',...
        'Tag','ipd',...
        'Position',[200,640,105,30],...
        'Callback',{@helpbutton_callback});
%     interpolate_text = uicontrol(mp,'Units','pixels','Style','text',...
%         'String','Numerical interpolation method','Position',[10,295,160,40]);
%     interpolate_pop = uicontrol(mp,'Style','popupmenu',...
%         'String',{'linear','spline','pchip','nearest'},...
%         'Position',[10,270,180,30],'Tag','interpop',...
%         'BackgroundColor','w',...
%         'Callback',{@interpop_callback});
%     intp_help = uicontrol(mp,'Units','pixels','Style','pushbutton',...
%         'String','Help',...
%         'Tag','intp',...
%         'Position',[200,274,105,30],...
%         'Callback',{@helpbutton_callback});
    im_panel = uipanel('Parent',mp,...
        'Title','Model system',...
        'Units','pixels',...
        'BackgroundColor','w',...
        'Position',[10,325,295,291]);
    im_ax = axes('Parent',im_panel,...
        'Units','pixels','Position',[10,10,275,251]);
    im=image(imread('model6.jpg'),'Parent',im_ax);
        axis off;
        %set(hObject,'Color','white');
    param_panel = uipanel('Parent',mp,...
        'Title','Loop parameters',...
        'Units','pixels',...
        'Position',[10,40,240,260]);
    param_table = uitable(param_panel,...
        'Data',dat,...
        'ColumnName',columnname,...
        'ColumnFormat',columnformat,...
        'ColumnEditable',columneditable,...
        'ColumnWidth',columnwidth,...
        'RowName',rowname_tab,...
        'Units','pixels',...
        'Position',[5,5,230,230],...
        'CellEditCallback',{@table_callback});
    ptab_help = uicontrol(mp,'Units','pixels','Style','pushbutton',...
        'String','Help',...
        'Tag','ptab',...
        'Position',[260, 45, 45, 247],...
        'Callback',{@helpbutton_callback});
    hold_check = uicontrol(f,'Style','checkbox',...
        'String','Keep the current plot',...
        'Value',0,...%'BackgroundColor','w',...
        'Position',[640,110,180,30],...
        'Callback',{@holdcheck_callback});
    simulate = uicontrol(f,'Units','pixels','Style','pushbutton',...
        'String','Simulate',...
        'Position',[675,70,90,30],...
        'Callback',{@simulatebutton_callback});
    mailauthor = uicontrol(f,'Units','pixels','Style','pushbutton',...
        'String','Contact Authors',...
        'Position',[665,30,110,30],...
        'Callback',{@mailauthor_callback});
    messages = uicontrol(f,'Units','pixels','Style','text',...
        'String',' ',...
        'Position',[570,10,315,20]);
    
    a1 = axes('Units','pixels','Tag','a1','Position',[90,90,450,405]);
    grid on
    a2 = axes('Units','pixels','Tag','a2','Position',[90,585,450,270]);
    grid on
    %a3 = axes('Position',[0.1,0.1,0.5,0.15]);
    %grid on
    
    
    % Initialize the code
    % Change units of the figure, so it resizes automatically
    set([f, bcgp, mp, ...
        inputdata_text, inputdata_pop, ipd_help, ...    
        param_panel, param_table, ptab_help, ...
        hold_check, simulate, mailauthor, messages, ...
        a1, a2, im_panel,im_ax],'Units','normalized');
    %    inputform_text,inputform_pop,ipf_help,...
    %    interpolate_text,interpolate_pop,intp_help,...
    
    % Generate initial data
    set([a1,a2],'FontSize',20)
    xlabel(a1,'Time t (ms)')
    ylabel(a1,'LGN response')
    ylabel(a2,'Cortical response')
    %ylabel(a3,'Response')
    
    % Make the GUI visible
    set(f,'Visible','on')
    
    function popup_callback(source, eventdata) 
        % Determine the selected input form
        str = get(source, 'String');
        val = get(source, 'Value');
        % Set current input to the selected input
        switch str{val};
            case 'Data file'
                in_struct.form = 'vec';
                set(inputdata_text,'String','Data file');
                set(inputdata_edit,'String','Provide a data file')
            case 'Function'
                in_struct.form = 'fnc';
                set(inputdata_text,'String','Function name/path');
                set(inputdata_edit,'String','Provide a function name')
            case 'In-line'
                in_struct.form = 'ilfnc';
                set(inputdata_text,'String','In-line function');
                set(inputdata_edit,'String','Provide an in-line function (use t for time)')
        end
        
    end
    
    function helpbutton_callback(source,eventdata)
        % Opens the respective help file
        str = get(source, 'Tag');
        hfile = ['helpFiles/' str '_help.pdf'];
        open(hfile);            
    end

    function datapop_callback(source, eventdata) 
        % Set data parameters
        str = get(source, 'String');
        val = get(source, 'Value');
        switch in_struct.form
            case 'vec'
                data=load(str{val});
                in_struct.t_in=data(:,1);
                in_struct.r_in=data(:,2);
                in_struct.data_file = str{val};
                clear data
                plot(a1,in_struct.t_in,in_struct.r_in, 'k--',...
                    'LineWidth',2);
                set(a1,'FontSize',20)
                legend(a1,'r_g')        
                ylabel(a1,'LGN response')
                [y1max,y1min] = find_max(a1);
                axis(a1,[0,max(in_struct.t_in),1.2*y1min,1.2*y1max])
                set([a1,a2],'XGrid','on','YGrid','on')
            case 'fnc'
                in_struct.h=user_string;
            case 'ilfnc'
        end
    end
    
    function interpop_callback(source, eventdata)
        % Determines the selected interpolation method
        str = get(source,'String');
        val = get(source,'Value');
        in_struct.interpol_method = str{val};
    end

    function table_callback(source, eventdata) %#ok<INUSL>
        % Update the parameter struct
        i=eventdata.Indices(1);
        val=str2double(eventdata.EditData);
        parname=cell2mat(rowname(i));
        lgn_struct.(parname)=val;
    end

    function holdcheck_callback(source, eventdata)
        val = get(source,'Value');
        if val == 0
            set([a1,a2],'NextPlot','replace')
            set(messages,'String',' ')
            ih1=1;
        elseif val ==1
            set([a1,a2],'NextPlot','add')
        end
    end
        
    function simulatebutton_callback(source, eventdata) %#ok<INUSD>
        % Simulate and plot
        
        % Create the waitbar
        w = waitbar(0,'Computing','Position',[300,300,288,60]);
        
        % load the first data file if there is no other data set 
        if strcmp(in_struct.data_file,' ')
            str = get(inputdata_pop,'String');
            val = get(inputdata_pop,'Value');
            defaultf=str{val};
            data=load(defaultf);
            in_struct.t_in=data(:,1);
            in_struct.r_in=data(:,2);
            in_struct.data_file = defaultf;
            clear data
        end
        
        % convolution of input with feedforward kernel
        [in_struct.t_g,in_struct.r_g]=filterMyInput(lgn_struct, in_struct);
        
        % evaluate the ode/dde governing the dynamics
        if (lgn_struct.Delta_fb==0)
            [t_rtmp,r_rtmp] = evaluate_ode(lgn_struct, in_struct);
        else
            [t_rtmp,r_rtmp] = evaluate_dde(lgn_struct, in_struct);
        end
        
        % Unpacking some parameters
        l_rON = lgn_struct.l_rON;
        l_cON = lgn_struct.l_cON;
        l_cOFF = lgn_struct.l_cOFF;
        
        % Setting the time vector
        t_r=in_struct.t_g;
        
        % smoothing the solution by interpolation        
        r_rn=interp1(t_rtmp,r_rtmp,t_r,in_struct.interpol_method);
        r_cON = (r_rn - l_cON) .* (r_rn - l_cON >= 0);
        r_cOFF = (-r_rn + l_cOFF) .* (-r_rn + l_cOFF >= 0);
        r_r = (r_rn - l_rON) .* (r_rn - l_rON >=0) + l_rON;
        clear t_rtmp r_rtmp r_cONtmp r_cOFFtmp
        % set the time vector for the cortical solution, no delay in
        % feedforward connection at the time...
        t_c = t_r;      % + lgn_struct.Delta_fb;
        % pack the solution and time vectors into the 'in_struct'
        in_struct.t_r=t_r;
        in_struct.t_c=t_c;
        in_struct.r_r=r_r;
        in_struct.r_cON=r_cON;
        in_struct.r_cOFF=r_cOFF;
        
        % dump parameters and in/out data to file tmp.mat
        save tmp.mat lgn_struct in_struct
        
        % unpack input vectors
        t_in=in_struct.t_in;
        r_in=in_struct.r_in;
            
        % ploting the results
        pl1 = plot(a1,[t_in;t_r(end)],[r_in;r_in(end)], '--',t_r,r_r,'-',...
            'LineWidth',2);
        set(a1,'FontSize',20)
        legend(a1,'r_g','r_r')        
        xlabel(a1,'Time t (ms)')
        ylabel(a1,'LGN response')
        [y1max,y1min] = find_max(a1);
        axis(a1,[0,max(t_r)-120,1.2*y1min,1.2*y1max])
        
        pl2 = plot(a2,[0; t_c],[0; r_cON],'-',...
            [0;t_c],[0; r_cOFF],'--',...
            'LineWidth',2);
        set(a2,'FontSize',20)
        legend(a2,'r_c^{ON}','r_c^{OFF}')
        ylabel(a2,'Cortical response')

        [y2max,y2min] = find_max(a2);
        if y2max == y2min
            axis(a2,[0,max(t_r),0,0.2])
        else
            axis(a2,[0,max(t_r)-120,1.2*y2min,1.2*y2max])
        end
        
%         plot(a3,[0; t_c],[0; r_cOFF],'k-',...
%             'LineWidth',2)
%         set(a3,'FontSize',20)
%         legend(a3,'r_c^{OFF}')
%         xlabel(a3,'Time t (ms)')
%         ylabel(a3,'Response')
%         [y3max,y3min] = find_max(a3);        
%         axis(a3,[0,max(t_r),1.2*y3min,1.2*y3max])
        
        set([a1,a2],'XGrid','on','YGrid','on')
        
        % setting colors for successive plots        
        %h1=get(a1,'Children');
        %ih1=length(h1)/2;
        linecol = ['k','b','r','c','m','y'];      
        if get(hold_check,'Value') == 1
            if ih1<length(linecol)
                ih1=ih1+1;
                col = linecol(ih1);                
            else
                set(messages,'String','All colors chosen, black will be used')
                col = 'k';
            end
        else 
            col = 'k';
        end
        set([pl1,pl2],'Color',col)        
        set(pl1(1),'Color','k')
        
        close(w)
        %set(messages,'String','Done');        
    end
    
    function mailauthor_callback(source, eventdata)
        web mailto:gaute.einevoll@nmbu.no
    end
end

function [yxmax,yxmin] =  find_max(ax)
    % function to find the max point of the lines in a plot
    hx=get(ax,'Children');
    yxraw = get(hx,'YData');
    if iscell(yxraw)
        m = size(yxraw);
        yxmax=0;
        yxmin=0;
        for icell = 1:m
            yxdata = cell2mat(yxraw(icell));
            yxmax_tmp = max(yxdata);
            yxmin_tmp = min(yxdata);
            if yxmax_tmp > yxmax
                yxmax=yxmax_tmp;
            end
            if yxmin_tmp < yxmin
                yxmin = yxmin_tmp;
            end
        end
    else
        yxmax=max(max(yxraw,[],2));
        yxmin=min(min(yxraw,[],2));
    end
end