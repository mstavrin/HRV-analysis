% 22 DEC 2009 -> EXTRACTING THE SEQUENCE OF FILES 

[filename1, pathname1]=uigetfiles('*.txt', 'Pick the dataset .mat to be analyzed!');
        if isequal(filename1, 0)|isequal(pathname1, 0)
            disp('User pressed cancel')
        else
            disp(['User selected ', fullfile(pathname1, filename1)])
        end
        cd(pathname1)
        a7=cell2struct(filename1, 'name', 1);
        for i=1:length(a7)
            current_name=a7(i,1).name;
            data=load(current_name)
     