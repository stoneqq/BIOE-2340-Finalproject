%Version - 1.10.2017
classdef NeuronGUI < handle
    properties
        nip
        handle
        editBoxes
        parameters
        nextActionTextbox
        ndA
        L
        nucleusBorder
        Cluster
        Single
        Small
        NClusterText
        PCluster
        cellBorder
        figaxis
        PCell
        NCellText
        ocbm
        ccnm
        cunm
        dirin
        dirout
        filein
        altparam
        h
        enableProcessOut
        enableProcessIn
        waitbar
        batch
        controlHandles
        flag
        parent
        batchInput
        batchOutput
        hSlider
        hControlpanel
        subPanel
        legend
        legend2
        badDataWarning = false; % flag indicating if the user has been notified about invalid parameter in editbox
        tableWindow
        pushButtonHandle
        
    end
    
    
    methods
        function handle=getHandle(ngui)
            if isempty(ngui.handle)%check if the figure window has be created
                ngui.handle=figure;
                handle=ngui.handle;
            else
                if ishandle(ngui.handle) %check if the figure window is currently open
                    handle=ngui.handle;
                else
                    ngui.handle=figure;
                    handle=ngui.handle;
                end
            end
            set(handle,  'name', 'Figure Window','numbertitle','off')
            set(handle, 'Units', 'normalized', 'Position', [ 0.32    0.05   0.636   0.8])
            ngui.subPanel = uipanel('Parent',handle,...
                'BackgroundColor','white',...
                'units', 'normalized',...
                'Position',[0.005,0.72,0.12,0.18]);
            ngui.legend = uicontrol('Parent', ngui.subPanel,...
                'Style','text',... %instruction for each parameter
                'units', 'normalized',...
                'FontSize', 10,...%11.5 - 7/6
                'BackgroundColor', [0.9, 0.9, 0.9],...
                'HorizontalAlignment','left',...
                'position',[0, 0.5, 1, 0.5]);
            ngui.legend2 = uicontrol('Parent', ngui.subPanel,...
                'Style','text',... %instruction for each parameter
                'units', 'normalized',...
                'FontSize', 10,...%11.5 - 7/6
                'BackgroundColor', [0.9, 0.9, 0.9],...
                'HorizontalAlignment','left',...
                'position',[0, 0, 1, 0.5]);
        end
        function ngui=NeuronGUI(varargin)
            %             addpath('..')
            p = mfilename('fullpath');
            indices = strfind(p,filesep);
            ngui.parent = p(1:indices(end)); %find GUI's parent dir(GAIN folder)
            addpath(ngui.parent) %add the GAIN directory into the path
            
            ngui.nip=NeuronImageProcessor; %create and store the image processor obj
            
            if nargin>0
                status = ngui.nip.readParametersFile(varargin{1});
                if ~isempty(status)
                    error(status)
                end
            end
            ngui.parameters=ngui.nip.getParameters;
            p = ngui.parameters(12);
            %             ngui.parameters.description
            [editBoxes, nextActionTextbox, controlpanelHandle, buttonHandles, sliderHandles,~,~,pushButtonHandle]=createControlPanel(ngui.parameters, ngui.nip.getActionName,@ngui.forwardButtonCallback, @ngui.backButtonCallback, @ngui.saveButtonCallback, @ngui.quitButtonCallback, @ngui.batchButtonCallback, @ngui.parameterButtonCallback, @ngui.sliderCallback, @ngui.editBoxCallback);
            ngui.editBoxes=editBoxes;
            ngui.hSlider = sliderHandles;
            ngui.nextActionTextbox=nextActionTextbox;
            ngui.handle=[];
            ngui.pushButtonHandle = pushButtonHandle;
            
            ngui.controlHandles = buttonHandles; %handles of buttons on control panel
            ngui.hControlpanel = controlpanelHandle; %handle of control panel window
            %hide the Figure Toolbar of the control panel window, name it,
            %and hide the number of the figure
            set(ngui.hControlpanel, 'menubar', 'none', 'name', 'Control Panel','numbertitle','off');
            
            %temp 10/28/16
            %show wait bar for the final image in interactive mode
            ngui.nip.showWaitBar(true)
            %/temp
        end
        function forwardButtonCallback(ngui,UIhandle,x)
            state=ngui.nip.getState();
            if state==NIPState.Ready
                [FileName, PathName]=uigetfile('*.*','Select the Nucleus Image File');
                nucleusImageFile=strcat(PathName,FileName);
                set(ngui.editBoxes(1),'string',nucleusImageFile);
                j=findjobj(ngui.editBoxes(1));%findjobj is a function from external resource
                j.setCaretPosition(length(get(ngui.editBoxes(1),'string')));
            end
            
            if ngui.badDataWarning%if user has been notified about invalid parameter in editbox
                ngui.badDataWarning=false;%flip the flag
            else
                for i=1:numel(ngui.editBoxes);
                    valueString=get(ngui.editBoxes(i),'string');
                    ngui.parameters(i).value=valueString;
                end
                status=ngui.nip.next(ngui.parameters);%call NIP
                updateUser(ngui,status);
            end
            %             cnt = 0;%count invalid values
            %             for i=1:numel(ngui.editBoxes)-1;
            %                 valueString=get(ngui.editBoxes(i+1),'string');
            %                 ngui.parameters(i+1).value=valueString;
            %                 number = str2num(ngui.parameters(i+1).value);
            %                 if isempty(number) || isnan(number) || numel(number) ~= 1
            %                     cnt = cnt+1;
            %                 end
            %
            %             end
            %             ngui.parameters(1).value = get(ngui.editBoxes(1),'string');%first edit box = image path
            %             if cnt == 0
            %             status=ngui.nip.next(ngui.parameters);
            %             updateUser(ngui,status);
            %             end
            %
            
            
        end
        
        function updateUser(ngui,status)
            state=ngui.nip.getState();
            %             fprintf('%s\n',char(state))
            if ~isempty(status)
                errordlg(status,'Error')
                return
            end
            %             fprintf('%s \n', char(state))
            
            %2/9/17
            if state == NIPState.Done
            %Disable the forward push button if the GUI is at the last
            %state
            set(ngui.pushButtonHandle, 'Enable','off')
            else
            set(ngui.pushButtonHandle, 'Enable','on')
            end
            
            switch(state)
                case NIPState.ReadImages
                    %                     I = ngui.nip.getCellImage();
                    %                     J = ngui.nip.getNucleusImage();
                    %                     h = figure(ngui.getHandle());
                    %                     imshow(J)
                    
                    %temp 10/27/16
                    rgb= createIntermediateImages(state,ngui.nip);
                    h = figure(ngui.getHandle());
                    imshow(rgb)
                    %/temp 10/27/16
                    %                     set(ngui.instructionTextbox, 'string', 'Original nucleus image')%instruction textbox on control panel
                    set(ngui.legend,'string',sprintf('%s\n%s\n%s','Original', 'Nucleus', 'Image'),'ForegroundColor', 'k')
                    
                    
                    %                     %Automatically update graphs when variables change
                    %                     for i = 2:numel(ngui.editBoxes)
                    %                         paraValue(i) = str2num(get(ngui.editBoxes(k+1),'String'));
                    %                     end
                    %                     linkdata on
                    
                case NIPState.SegmentedNucleusImageOnce
                    %                     I=ngui.nip.getFirstNucleusMask();
                    %                     J=ngui.nip.getNucleusImage();
                    %                      rgb = addBorder(J, I, [0, 0, 1]);
                    %                     figure(ngui.getHandle());
                    %                     imshow(rgb)
                    
                    %temp 10/27/16
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle());
                    imshow(rgb)
                    %/temp 10/27/16
                    
                    figure(ngui.getHandle());%without this command, the subPanel will not be shown (?)
                    %                    set(ngui.instructionTextbox, 'string', sprintf('blue - nuclei')) %instruction textbox on control panel
                    set(ngui.legend, 'String', sprintf('%s\n%s', 'Blue:', 'nuclei'), 'ForegroundColor', 'b')
                case NIPState.SegmentedNucleusImageTwice
                    %                     I=ngui.nip.getSecondNucleusMask();
                    %                     J=ngui.nip.getNucleusImage();
                    %                     figure(ngui.getHandle);
                    %                     imshow(rgb)
                    %temp 10/27/16
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle());
                    imshow(rgb)
                    %/temp 10/27/16
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s', 'Blue:', 'nuclei'), 'ForegroundColor', 'b')
                case NIPState.OpenedNucleusMask
                    %                     I=ngui.nip.getOpenedNucleusMask;
                    %                     J=ngui.nip.getNucleusImage();
                    %                     rgb = addBorder(J, I, [0, 0, 1]);
                    %                     figure(ngui.getHandle);
                    %                     imshow(rgb)
                    %temp 10/27/16
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle());
                    imshow(rgb)
                    %/temp 10/27/16
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s', 'Blue:', 'nuclei'), 'ForegroundColor', 'b')
                case NIPState.IdentifiedNucleusClusters
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    set(ngui.legend, 'String', sprintf('%s\n%s', 'Cyan:', 'nuclei clusters'), 'ForegroundColor', [0, 0.8, 0.8])
                    %If no cluster being identified, notify the users
                    %(2/9/17)
                    
                case NIPState.CalculatedNominalMeanNucleusArea
                    [rgb, cluster] = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s', 'Cyan:', 'nuclei clusters'), 'ForegroundColor', [0, 0.8, 0.8])
                    %legend2: if there are clusters, show nothing; if no
                    %cluster, report.
                    set(ngui.legend2, 'String', sprintf('%s', cluster), 'ForegroundColor', [0, 0, 0])
                case NIPState.CalculatedMinNucleusArea
                    [rgb, cluster] = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s','Magenta:', sprintf('%s\n%s', 'Nuclei too small', 'to be accepted')), 'ForegroundColor', 'magenta')
                    %legend2: if there are clusters, show nothing; if no
                    %cluster, report.
                    set(ngui.legend2, 'String', sprintf('%s', cluster), 'ForegroundColor', [0, 0, 0])
                case NIPState.SegmentedCells
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    set(ngui.legend, 'String', sprintf('%s\n%s\n%s','Cell', 'Image', 'Opened'))
                    set(ngui.legend2, 'String', 'Red: cell bodies and neurites', 'ForegroundColor', 'r')
                case NIPState.SeparatedBodiesFromNeurites
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s','Red:', 'cell bodies'), 'ForegroundColor', 'r')
                    set(ngui.legend2, 'String', sprintf('%s\n%s','Green:', 'neurites'), 'ForegroundColor', [0, 0.9, 0])
                case NIPState.ResegmentedNeurites
                    rgb = createIntermediateImages(state,ngui.nip);
                    figure(ngui.getHandle);
                    imshow(rgb)
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s','Green:', 'connected neurites'), 'ForegroundColor', [0, 0.9, 0])
                    set(ngui.legend2, 'String', sprintf('%s\n%s', 'Yellow:', 'unconnected neurites'), 'ForegroundColor', [0.78, 0.78, 0])
                case NIPState.ResegmentedNeuriteEdges  %3rd Neurite Segmentation - added on 6/17/16
                    %                     cnm = ngui.nip.getThirdConnectedNeuriteMask();
                    %                     unm = ngui.nip.getThirdUnconnectedNeuriteMask();
                    %                     I = ngui.nip.getCellImage();
                    %                     rgb = addBorder(I, ngui.ocbm, [1, 0, 0]);
                    %                     rgb = addBorder(rgb, cnm, [0, 1, 0]);
                    %                     rgb = addBorder(rgb, unm, [1, 1, 0]);
                    %                     rgb = insertText(rgb,[ngui.PCell(2,:)',ngui.PCell(1,:)'],ngui.NCellText,'TextColor','white','FontSize',28, 'BoxOpacity',0);
                    
                    %temp 10/27/16
                    rgb = createIntermediateImages(state,ngui.nip);
                    %/temp
                    figure(ngui.getHandle);
                    imshow(rgb)
                    
                    set(ngui.legend, 'String', sprintf('%s\n%s','Green:', 'connected neurites'), 'ForegroundColor', [0, 0.9, 0])
                    set(ngui.legend2, 'String', sprintf('%s\n%s', 'Yellow:', 'unconnected neurites'), 'ForegroundColor', [0.78, 0.78, 0])
                case NIPState.ClosedNeuriteMask
                    %                     ngui.ccnm = ngui.nip.getClosedConnectedNeuriteMask();
                    %                     ngui.cunm = ngui.nip.getClosedUnconnectedNeuriteMask();
                    %                     I = ngui.nip.getCellImage();
                    %                     rgb = addBorder(I, ngui.ocbm, [1, 0, 0]);
                    %                     rgb = addBorder(rgb, ngui.ccnm, [0, 1, 0]);
                    %                     rgb = addBorder(rgb, ngui.cunm, [1, 1, 0]);
                    %                     rgb = insertText(rgb,[ngui.PCell(2,:)',ngui.PCell(1,:)'],ngui.NCellText,'TextColor','white','FontSize',28, 'BoxOpacity',0);
                    
                    %temp 10/27/16
                    rgb = createIntermediateImages(state,ngui.nip);
                    %/temp
                    
                    figure(ngui.getHandle);
                    imshow(rgb)
                    set(ngui.legend, 'String', sprintf('%s\n%s','Green:', 'connected neurites'), 'ForegroundColor', [0, 0.9, 0])
                    set(ngui.legend2, 'String', sprintf('%s\n%s', 'Yellow:', 'unconnected neurites'), 'ForegroundColor', [0.78, 0.78, 0])
                    
                    
                    
                case NIPState.Done
                    %display the final image
                    rgb = createIntermediateImages(state,ngui.nip);
                
                    figure(ngui.getHandle);
                    imshow(rgb)
                    set(ngui.legend, 'String', sprintf('%s\n%s','Green:', 'Long neurites'), 'ForegroundColor', [0, 0.9, 0])
                    set(ngui.legend2, 'String', sprintf('%s\n%s', 'Blue:', 'Short neurites'), 'ForegroundColor', 'b')
                    
                    %temp 11/1/16
                    %create a data table and display the table on a figure
                    %window
                    [dHeader,dataCell] = createDataTable(ngui.nip);
                    fTable = figure;
                    tHandle = uitable(fTable);
                    set(tHandle, 'Data', dataCell, 'Units', 'normalized','Position' , [0,0,1,1], 'RowName', ({}), 'ColumnName',dHeader);
                   %2/9/17
                   msgbox('Image processing is complete.')
                    
                    
                    
                otherwise
                    error('[NeuronGUI.updateUser] Unexpected State: %s', char(state));
            end
            
            %check if ngui.nip.getActionName will go beyond one line
            
            %edited on 7/5
            if length(ngui.nip.getActionName) < 26 %Hard-coded. Don't know how to check if a line in textbox is full
                actionStr = sprintf('\n%s',ngui.nip.getActionName); %if only one line, move it to the bottom of the textbox
            else
                actionStr = ngui.nip.getActionName;
            end
            set(ngui.nextActionTextbox,'string',actionStr);
            
            %             set(ngui.instructionTextbox,'string',ngui.nip.getActionName);%instruction of each parameter;
            ngui.parameters=ngui.nip.getParameters;
            for i=1:numel(ngui.parameters)
                if ngui.parameters(i).active
                    enbl='on';
                    %                     ngui.parameters(i).description % descriptions are empty now
                else
                    enbl='off';
                end
                set(ngui.editBoxes(i),'Enable',enbl);
                if i>1
                    set(ngui.hSlider(i-1),'Enable',enbl);
                end
            end
            
        end
        
        function backButtonCallback(ngui,UIhandle,x)
            
            if ngui.badDataWarning
                ngui.badDataWarning = true;
                
            else
                %                 ngui.badDataWarning
                %                 %If any parameter value is empty, cannot go back
                %             for i=2:numel(ngui.editBoxes)
                %                 if isempty(get(ngui.editBoxes(i),'String'))
                %                     warndlg('Empty parameter value is not allowed. Please select a value')
                %                     return
                %                 end
                %             end
                for i=1:numel(ngui.editBoxes);
                    valueString=get(ngui.editBoxes(i),'string');
                    ngui.parameters(i).value=valueString;
                end
                status=ngui.nip.back(ngui.parameters);
                updateUser(ngui,status);
                
            end
        end
        
        
        function quitButtonCallback(ngui,UIhandle, x)
            %             close all %need to check if that is sufficient
            if ishandle(ngui.handle)
                close(ngui.handle)%close the figure window
            end
            close(ngui.hControlpanel)%close the control panel window
            
            rmpath(ngui.parent)   %remove the GAIN directory when the user quits the program
        end
        function batchButtonCallback(ngui,UIhandle, x)% "Exit to batch processing" button
            %temp 10/28/16
            %Hide wait bar for the final image in batch mode
            ngui.nip.showWaitBar(false)
            %/temp
            ngui.batch = figure;
            %hide the figure toolbar
            set(ngui.batch, 'menubar', 'none', 'name', 'Settings for Batch Processing','numbertitle','off');
            ngui.batch,%check if it is necessary
            %temp
            ngui.batch.Position(4) = 5/6*ngui.batch.Position(4); %decrease the height of the batch window - 9.18
            wBatch = 560; %the wdith of batch window - calibrated
            hBatch = 350; %the height of batch window - calibrated
            
            bottom = 0;
            processButtonHandle = uicontrol('Style', 'pushbutton', 'String', 'Process',...
                'units', 'normalized',...
                'Position', [410/wBatch bottom 150/wBatch 50/hBatch], 'Callback', @ngui.processButtonCallback, ...
                'Enable', 'Off');
            uicontrol('Style', 'pushbutton', 'String', 'Cancel',...%original name: Back to Control Panel - 9.16.2016
                'units', 'normalized',...
                'Position', [0 bottom 150/wBatch 50/hBatch], 'Callback', {@ngui.controlPanelButtonCallback});%close batch wondow and back to controlpanel
            bottom_out = bottom + 125/hBatch;
            uicontrol('Style', 'pushbutton', 'String', 'Output Directory',...
                'units', 'normalized',...
                'Position', [0 bottom_out 100/wBatch 25/hBatch], 'Callback', {@ngui.outputButtonCallback, processButtonHandle});
            bottom_infile = bottom +215/hBatch;
            bottom_indir = bottom_infile+75/hBatch;
            uicontrol('Style', 'pushbutton', 'String', 'Input Files',...
                'units', 'normalized',...
                'Position', [0 bottom_infile 100/wBatch 25/hBatch], 'Callback', {@ngui.inputFileButtonCallback, processButtonHandle});
            uicontrol('Style', 'pushbutton', 'String', 'Input Directory',...
                'units', 'normalized',...
                'Position', [0 bottom_indir 100/wBatch 25/hBatch], 'Callback', {@ngui.inputButtonCallback, processButtonHandle});
            
            %             set(ngui.batch, 'Position', [0.01 0.06 wControlPanel_rel hgtControlPanel_rel]);%Auto resize the window to fit the user's screen
            uicontrol('Style','text',... % a sign between "input dir" and "input files" buttons
                'units', 'normalized',...
                'string','Or',...
                'FontUnits', 'normalized',...
                'FontSize', 0.8,...
                'position', [35/wBatch bottom_infile+40/hBatch 30/wBatch 20/hBatch])
            
            ngui.flag = 0; %set the flag to be 0
            
            ngui.batchOutput = uicontrol('Style', 'edit', 'String', ' ',  'units', 'normalized',...
                'Position', [110/wBatch bottom_out-25/hBatch 400/wBatch 50/hBatch],'HorizontalAlignment','left', ...
                'Max', 100, 'Enable', 'on');
            jOutput=findjobj(ngui.batchOutput,'nomenu'); %get the UIScrollPane container
            jOutput=jOutput.getComponent(0).getComponent(0);
            set(jOutput,'Editable',0);
            
            ngui.batchInput = uicontrol('Style', 'edit',  'units', 'normalized',...
                'string', ' ', ...
                'Max', 100,...
                'Position', [110/wBatch bottom_infile 400/wBatch 100/hBatch],'HorizontalAlignment','left', ...
                'Enable', 'on');%off 6/28
            
            jInput=findjobj(ngui.batchInput,'nomenu'); %get the UIScrollPane container
            jInput=jInput.getComponent(0).getComponent(0);
            set(jInput,'Editable',0);
            
            
            for k = 1:length(ngui.controlHandles)
                set(ngui.controlHandles{k},'Enable','off')
            end
            
            %figure close call back
            % close batch processing window when a user clicks "back to
            % control panel"
            set(ngui.batch, 'CloseRequestFcn',@ngui.controlPanelButtonCallback)
            
        end
        function controlPanelButtonCallback(ngui,UIhandle, x) %back to control panel(cancel button)
            %                 error('Terminate batch processing.')
            %temp 10/28/16
            %show wait bar for the final image in interactive mode
            ngui.nip.showWaitBar(true)
            %/temp
            ngui.flag = 1;%clicking the button changes the flag - which terminates the processing
            if ishandle(ngui.batch)
                close(ngui.batch)   %close batch processing window
            end
            %             if ishandle(ngui.waitbar), close(ngui.waitbar), end %close the waitbar window if it is open
            %             close(ngui.waitbar)
            for k = 1:length(ngui.controlHandles)
                set(ngui.controlHandles{k},'Enable','on')
            end
        end
        
        function outputButtonCallback(ngui, UIhandle, x, processButtonHandle)
            ngui.dirout=uigetdir('*.*','Store Data');
            
            if ischar(ngui.dirout) %if the directory is valid (is a string)
                set(ngui.batchOutput,'string', ngui.dirout);
            else
                set(ngui.batchOutput,'string', '');%if not, don't show it on the editbox
            end
            
            
            if ~isempty(ngui.dirout) && ischar(ngui.dirout) %matlab returns numerical 0 when nothing is selected
                ngui.enableProcessOut = 'On';
            else
                ngui.enableProcessOut = 'Off';
            end
            %             if ngui.dirout == 0, ngui.enableProcessOut = 'Off'; end %Click output  button but not select anything
            %             ngui.enableProcessOut = 'On'; %do we need to check if ~isempty(ngui.dirin)??
            if strcmp(ngui.enableProcessOut, 'On') && strcmp(ngui.enableProcessIn, 'On')
                enableProcess = 'On';
            else
                enableProcess = 'Off';
            end
            set(processButtonHandle,'Enable',enableProcess);
        end
        
        function inputButtonCallback(ngui, UIhandle, x, processButtonHandle)%get directory input
            %reset ngui.filein, since ngui.dirin is going to be used, not ngui.filein
            ngui.filein = [];
            ngui.dirin=uigetdir('*.*','Image File');
            if ischar(ngui.dirin)%if the directory is valid (is a string)
                set(ngui.batchInput, 'string', ngui.dirin)
            else set(ngui.batchInput,'string', '');%if not, don't show it on the editbox
            end
            
            if ~isempty(ngui.dirin) && ischar(ngui.dirin) %matlab returns numerical 0 when nothing is selected
                ngui.enableProcessIn = 'On';
            else
                ngui.enableProcessIn = 'Off';
            end
            if strcmp(ngui.enableProcessOut, 'On') && strcmp(ngui.enableProcessIn, 'On')
                enableProcess = 'On';
            else
                enableProcess = 'Off';
            end
            set(processButtonHandle,'Enable',enableProcess);
        end
        
        function inputFileButtonCallback(ngui, UIhandle, x, processButtonHandle)%get files input
            %reset ngui.dirin, since ngui.filein is going to be used, not ngui.dirin
            ngui.dirin = [];
            [file,path]=uigetfile('*.*','Image File', 'MultiSelect','on');%same path for the files
            if iscell(file) %if more then 1 file is selected (cell array)
                ngui.filein = cell(1, length(file));
                for i = 1: length(file)
                    ngui.filein{i} = strcat(path,file{i});
                end
                fileNames = sprintf('%s\n',ngui.filein{:});%convert a cell of char vectors to a multi-line char vector
            else       %if only one file is selected (a character vector)
                ngui.filein = strcat(path,file);
                fileNames = ngui.filein;
            end
            set(ngui.batchInput, 'string', fileNames)
            
            if ~isempty(ngui.filein)
                ngui.enableProcessIn = 'On';
            else
                ngui.enableProcessIn = 'Off';
            end
            if strcmp(ngui.enableProcessOut, 'On') && strcmp(ngui.enableProcessIn, 'On')
                enableProcess = 'On';
            else
                enableProcess = 'Off';
            end
            set(processButtonHandle,'Enable',enableProcess);
        end
        
        
        function parameterButtonCallback(ngui, UIhandle, x) %call back for button "Open Parameter File"
            [FileName, PathName] = uigetfile('*.*', 'Select Parameters File');
            ngui.altparam = ngui.nip.readParametersFile(strcat(PathName, FileName));%check the status
            if ~isempty(ngui.altparam)
                errordlg(ngui.altparam,'Error')
                return
            end
            
            ngui.parameters=ngui.nip.getParameters;
            
            [editBoxesnew, hSlidernew] = updateControlPanel(ngui.editBoxes, ngui.parameters, ngui.h, ngui.hSlider);
            ngui.editBoxes=editBoxesnew;
            ngui.hSlider=hSlidernew;
        end
        
        function processButtonCallback(ngui, UIhandle, x)
            if ~isempty(ngui.dirin)  %if the input is a directory
                list=dir(ngui.dirin); %get the list of the file names under the directory
                namelist = cell(1,(length(list)-2));
                for i = 1:length(namelist)
                    namelist{i} = list(i+2).name;
                    namelist{i} = strcat(ngui.dirin, filesep, namelist{i});
                end
            elseif ~isempty(ngui.filein) %if the input is files
                namelist = ngui.filein;
            else
                error('The input cannot be empty')
            end
            %modified on 1/10/17
            tElapsed = 216;%initial guess of the processing time for an image = 216s (based on tests)
            
            if iscell(namelist) % if the input is more than one file
                currentWait = round(tElapsed*(length(namelist))/60,1);
                %2/9/17
%                 ngui.waitbar = waitbar(1/length(namelist), ['Processing Image 1 of ' num2str(length(namelist)) ' Approximate Time: ' num2str(currentWait) 'minutes'])
                ngui.waitbar = waitbar(0, ['Processing Image 1 of ' num2str(length(namelist)) ' Approximate Time: ' num2str(currentWait) 'minutes'])
                for i = 1:length(namelist)
                    if ngui.flag == 1, break,end %check the flag. If button "back to control panel" is clicked, flag = 1, otherwise flag = 0
                    tStart = tic;
                    ngui.nip.oneProcess(namelist{i}, ngui.dirout);
                    tElapsed = toc(tStart);
                    currentWait = round(tElapsed*(length(namelist)-i)/60,1);
                    
                    if i<length(namelist)
                        %indicate the ith image has already been processed
                        waitbar((i)/length(namelist), ngui.waitbar, ['Processing Image ' num2str(i+1) ' of ' num2str(length(namelist)) ' Approximate Time: ' num2str(currentWait) 'minutes'])
                    else
                        close(ngui.waitbar)
                    end
                end
            else  %if the input is only one file
                currentWait = round(tElapsed/60,1);
                ngui.waitbar = waitbar(0, ['Processing Image 1 of 1  Approximate Time: ' num2str(currentWait) 'minutes'])
                tStart = tic;
                ngui.nip.oneProcess(namelist, ngui.dirout);
                tElapsed = toc(tStart);
                close(ngui.waitbar)
            end
            
            %2/9/17
            %When the image processing task is complete, create a message
            %box to notify the user.
            msgbox('Batch image processing is complete.');
        end
        
        function saveButtonCallback(ngui, UIhandle, x)
            
            [FileName, PathName]=uiputfile('*.*','Save parameters as');
            %Do nothing if the file name or path name is a character array
            %when the user cancels the pop-up window, FileName and PathName
            %will be 0 (double).
            if ischar(FileName) && ischar(PathName)
                parameterData=strcat(PathName, FileName);
                %get status
                
                status = ngui.nip.writeParametersFile(parameterData);
                if ~isempty(status)
                    errordlg(status,'Error')
                    return
                end
                
            end
        end
        
        function sliderCallback(ngui, UIhandle, x)
            %(1) Update edit boxes
            num = length(ngui.hSlider);%num of sliders
            sliderValue = cell(1,num);
            for k = 1:num
                sliderValue{k} = num2str(get(ngui.hSlider(k),'Value'));
                set(ngui.editBoxes(k+1),'String', sliderValue{k})
            end
            
            
            
            
        end
        
        function editBoxCallback(ngui, UIhandle, x, paraName,subtype, editBoxID)
            newValueStr = get(UIhandle, 'String');
            newValueNum = str2num(newValueStr);
            status = subtype.check(newValueNum);
            if isempty(status)
                ngui.badDataWarning = false;%no bad data
                textValue = newValueNum;
                sliderID = editBoxID-1;
                set(ngui.hSlider(sliderID),'Value', textValue)
                if (newValueNum > ngui.hSlider(sliderID).Max) && (newValueNum ~= inf)
                    set(ngui.hSlider(sliderID),'Max', newValueNum)
                end
            else
                ngui.badDataWarning = true;%badData has been detected
                errordlg(strcat('Parameter',{' '} , paraName, {': '}, status))
            end
            %             if isempty(newValueNum)
            %                 errordlg(strcat('The value of parameter', paraName, 'is invalid. It must be a number or Inf.'))
            %                 ngui.badDataWarning = true;%badData has been detected
            %             elseif  isnan(newValueNum)
            %                 errordlg(strcat('The value of parameter', paraName, ' is Not a Number. It must be a number or Inf.'))
            %                 ngui.badDataWarning = true;%badData has been detected
            %             elseif numel(newValueNum) ~= 1
            %                 errordlg(strcat('The value of parameter', paraName, ' is not a single number'))
            %                 ngui.badDataWarning = true;%badData has been detected
            %             else
            %                 ngui.badDataWarning = false;%no bad data
            %                 textValue = newValueNum;
            %                 set(ngui.hSlider,'Value', textValue)
            %                 if (newValueNum > ngui.hSlider.Max) && (newValueNum ~= inf)
            %                     set(ngui.hSlider,'Max', newValueNum)
            %                 end
            %                 ngui.badDataWarning
            %             num = length(ngui.hSlider);%num of sliders
            %             textValue = zeros(1,num);
            %             for k = 1:num
            %                 newValueStr = get(ngui.editBoxes(k+1),'String'); %new parameter value string
            %                 newValueNum = str2num(newValueStr); %new parameter value number
            %                 if isempty(newValueNum)
            %                     errordlg(strcat('The value of parameter', num2str(k), 'is invalid. It must be a number or Inf.'))
            %                     ngui.badDataWarning = true;%badData has been detected
            %                 elseif  isnan(newValueNum)
            %                     errordlg(strcat('The value of parameter', num2str(k), ' is Not a Number. It must be a number or Inf.'))
            %                     ngui.badDataWarning = true;%badData has been detected
            %                 elseif numel(newValueNum) ~= 1
            %                     errordlg(strcat('The value of parameter', num2str(k), ' is not a single number'))
            %                     ngui.badDataWarning = true;%badData has been detected
            %                 else
            %                     ngui.badDataWarning = false;%no bad data
            %                     textValue(k) = newValueNum;
            %                     set(ngui.hSlider(k),'Value', textValue(k))
            %                     if (newValueNum > ngui.hSlider(k).Max) && (newValueNum ~= inf)
            %                     set(ngui.hSlider(k),'Max', newValueNum)
            %                     end
            %                  end
            
            
        end
        
        
        
    end
    
    
    %         function
    %         end
end








