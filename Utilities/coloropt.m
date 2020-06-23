%% color functions
function options = coloropt(cmod,linesty,marker,alpha)

    options.line_width = 2;
    options.error      = 'std';
    
    options.alpha      = alpha;
    options.line_style = linesty;
    options.marker     = marker;

    switch cmod
        case 1
            options.color_area = [128 193 219]./255;    % Blue theme
            options.color_line = [ 52 148 186]./255;
        case 2
            options.color_area = [243 169 114]./255;    % Orange theme
            options.color_line = [236 112  22]./255;
        case 3
            options.color_area = [ 31 241 137]./255;    % Green theme
            options.color_line = [  1 196   1]./255;
    end

end