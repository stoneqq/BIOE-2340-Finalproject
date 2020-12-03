%% probing the functions and variables in GAIN
% resegmentNeurites

ngui.nip=NeuronImageProcessor; %create and store the image processor obj
ngui.dirout="E:\Stoneqq\Workspace\Education\Pitt\Course\BIOE2340-Intro to Medical Imaging and Image Analysis\Project\GAIN\Test";
ngui.parameters=ngui.nip.getParameters;
state=ngui.nip.getState();

if state==NIPState.Ready
    [FileName, PathName]=uigetfile('*.*','Select the Nucleus Image File');
    nucleusImageFile=strcat(PathName,FileName);
    ngui.nip.oneProcess_nooutput(nucleusImageFile, ngui.dirout);
    ngui.nip.readImageFile();
    rgb= ngui.nip.getCellImage();
    imshow(rgb)
    rgb2= ngui.nip.getNucleusImage();
    imshow(rgb2)
end
ngui.nip.segmentCellBodies();
ngui.nip.getCellThresh1();
[rgb,cluster] = createIntermediateImages(NIPState.SegmentedCells,ngui.nip);
%        ngui.nip.isolateCellBodies();
%        ngui.nip.resegmentNeurites();
%        ngui.nip.resegmentNeuriteEdges();
%        ngui.nip.secondNeuriteResegmentation();
       
