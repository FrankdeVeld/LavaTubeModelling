%% Function for writing a .model file, needed foor IGMAS modelling
% Function to write a basic .model file, according to the previously
% defined model limits
function Write_basic_model_file(Rock_density, Cavity_density, Model_name, Model_limits, Cavity_limits, Num_intermediate_sections) 
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work
    File_name = [Model_name,'_Basic_model.model'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % Result in proper encoding
    fprintf(fid, '<geodata name="Test Model">\n'); % Open the geodata class
    
    % Define a reference body
    fprintf(fid, '<property name="body" value="reference"> \n    <property name="density" units="g/cm³" value="0.0"></property> \n <color red="0.5019608" green="0.5019608" blue="0.5019608"></color> \n</property> \n' );

    % Define another, rocky body. Note: density can be changed!
    String_rock = ['<property name="body" value="stone"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="0.6601" green="0.6601" blue="0.6601"></color> \n</property> \n'];
    String_rock = sprintf(String_rock,Rock_density);
    fprintf(fid, String_rock);

    % Define the cavity body. Density should be 0, but can be changed
    String_cavity = ['<property name="body" value="cavity"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="0" green="1" blue="0"></color> \n</property> \n'];
    String_cavity = sprintf(String_cavity,Cavity_density);
    fprintf(fid, String_cavity);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Section 0
    
    Section_text = ['<geometry>' ...
      '<cross_section name="0" '...
      'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
    Section_text = sprintf(Section_text,Model_limits(1),Model_limits(3),Model_limits(2),Model_limits(3));

    Temp_text = ['<vertex id="0" x="%4.2f" z="%4.2f"></vertex> \n' ...
    '<vertex id="1" x="%4.2f" z="%4.2f"></vertex> \n' ....
    '<vertex id="2" x="%4.2f" z="%4.2f"></vertex> \n'...
    '<vertex id="3" x="%4.2f" z="%4.2f"></vertex> \n'];
    Temp_text = sprintf(Temp_text,Model_limits(2)-Model_limits(1),Model_limits(5),0,Model_limits(5),0,Model_limits(6),Model_limits(2)-Model_limits(1),Model_limits(6));
    Section_text = [Section_text,Temp_text];

    Section_text = [Section_text,'<entry type="polygon" id_list="0 1 2 3 "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n']; 

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Intermediate Sections 
    index = 3;
    
    for a=1:Num_intermediate_sections
        Section_text = ['<geometry>' ...
          '<cross_section name="%d" '...
          'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
        Section_text = sprintf(Section_text,a,Model_limits(1),a*(Cavity_limits(4)-Cavity_limits(3))/(Num_intermediate_sections + 1)+ Cavity_limits(3),Model_limits(2),a*(Cavity_limits(4)-Cavity_limits(3))/(Num_intermediate_sections + 1)+ Cavity_limits(3));

        Temp_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ...
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ....
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'...
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
        Temp_text = sprintf(Temp_text,index + 1, Model_limits(2)-Model_limits(1),Model_limits(5),index + 2,0,Model_limits(5),index + 3,0,Model_limits(6),index + 4,Model_limits(2)-Model_limits(1),Model_limits(6));
        Section_text = [Section_text,Temp_text];
        
        Temp_text = ['<entry type="polygon" id_list="%d %d %d %d "> \n' ...
             '<property name="body" value="stone"></property> \n'...
          '</entry> \n'];
        Temp_text = sprintf(Temp_text,index+1,index+2,index+3,index+4);
        Section_text = [Section_text,Temp_text];
        
        index = index + 4;
        
        Section_text = [Section_text,'</geometry> \n'];
        fprintf(fid,Section_text);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Last section
    
    Section_text = ['<geometry>' ...
      '<cross_section name="%d" '...
      'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
    Section_text = sprintf(Section_text,a+1,Model_limits(1),Model_limits(4),Model_limits(2),Model_limits(4));

    Temp_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ...
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ....
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'...
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
    Temp_text = sprintf(Temp_text,index + 1, Model_limits(2)-Model_limits(1),Model_limits(5),index + 2, 0,Model_limits(5),index + 3,0,Model_limits(6),index + 4,Model_limits(2)-Model_limits(1),Model_limits(6));
    Section_text = [Section_text,Temp_text];
    
    Temp_text = ['<entry type="polygon" id_list="%d %d %d %d "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n'];
    Temp_text = sprintf(Temp_text,index+1,index+2,index+3,index+4);
    Section_text = [Section_text, Temp_text];

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);


    fprintf(fid, '</geodata>\n' );  % Close the geodata class, and with this the file
    fclose(fid);
end
