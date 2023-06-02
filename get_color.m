function color = get_color(name)
% Description 

switch name
    case 'bondiBlue1' % from https://colorswall.com/palette/25326
        color = [230, 244, 248]/255; % lightest bondi blue 
        
    case 'bondiBlue2'
        color = [204, 234, 240]/255;
        
    case 'bondiBlue3'
        color = [179, 223, 233]/255;
        
    case 'bondiBlue4'
        color = [153, 213, 226]/255;
        
    case 'bondiBlue5'
        color = [128, 202, 219]/255;
        
    case 'bondiBlue6'
        color = [102, 191, 211]/255;
        
    case 'bondiBlue7'
        color = [77, 181, 204]/255;
        
    case 'bondiBlue8'
        color = [51, 170, 197]/255;
        
    case 'bondiBlue9'
        color = [25, 160, 189]/255;
        
    case 'bondiBlue10'
        color = [0 149 182]/255;
        
    case 'bondiBlue11'
        color = [0, 134, 164]/255;
        
    case 'bondiBlue12'
        color =  [0, 119, 146]/255;
        
    case 'bondiBlue13'
        color = [0, 104, 127]/255;
        
    case 'bondiBlue14'
        color = [0, 89, 109]/255;
        
    case 'bondiBlue15'
        color = [0, 75, 91]/255;
        
    case 'bondiBlue16'
        color =  [0, 60, 73]/255; % darkest bondi blue
        
    case 'darkGrey'
        color = [.4 .4 .4];
        
    case 'salmon1' % from https://www.color-hex.com/color-palette/2539
        color = [255,186,186]/255;
        
    case 'salmon2'
        color = [255,123,123]/255;
        
    case 'salmon3'
        color = [255,82,82]/255;
        
    case 'red'
        color = [255,0,0]/255;
        
    case 'darkRed'
        color = [167,0,0]/255;
        
    case 'matcha1'
        color = [172, 209, 175]/255; % from https://colorswall.com/palette/45134
        
    case 'matcha2'
        color = [177, 221, 158]/255;
        
    case 'matcha3'
        color = [151, 207, 138]/255;
        
    case 'matcha4'
        color = [122, 159, 121]/255;
        
    case 'matcha5'
        color = [49, 90, 57]/255;
        
    case 'matcha6'
        color = [49, 94, 38]/255;
        
    case 'orange1'
        color = [255, 175, 122]/255; % from https://colorswall.com/palette/1322/
        
    case 'orange2'
        color = [255, 157, 92]/255;
        
    case 'orange3'
        color = [255, 139, 61]/255;
        
    case 'orange4'
        color = [255, 120, 31]/255;
        
    case 'orange5'
        color = [255, 102, 0]/255;
        
    case 'gold1'
        color = [229, 207, 135]/255; % from https://colorswall.com/palette/1322/
        
    case 'gold2'
        color = [221, 191, 95]/255;
        
    case 'gold3'
        color = [216, 183, 75]/255;
        
    case 'gold4'
        color = [212, 175, 55]/255;
        
    case 'gold5'
        color = [170, 140, 44]/255;
        
    case 'gold6'
        color = [127, 105, 33]/255;
end

end