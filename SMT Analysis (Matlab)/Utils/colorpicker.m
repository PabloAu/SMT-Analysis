function varargout = colorpicker(varargin)
% colorpicker gets RGB values from mouse clicks anywhere on a figure. Press
% Return to exit. 
% 
%% Syntax
% 
%  colorpicker
%  colorpicker(256)
%  colorpicker(...,TextProperty,TextValue) 
%  RGB = colorpicker(...)
%  [R,G,B] = colorpicker(...) 
% 
%% Description 
%
% colorpicker gets RGB values of anywhere you click on the current figure. 
% Values are temporarily printed at the location of the click until you 
% click somewhere else or press Return to exit. 
%
% colorpicker(256) scales values from 0 to 256 rather than the Matlab 
% default method of scaling from 0 to 1. 
%
% colorpicker(...,TextProperty,TextValue) specifies text formatting for the
% temporary RGB labels. 
% 
% RGB = colorpicker(...) returns RGB values of all clicks while the color 
% picker is open.  If only one output is requested, RGB values a an N-by-3
% array, where N is the number of clicks. 
% 
% [R,G,B] = colorpicker(...) returns RGB values as three individual N-by-1 
% arrays.  I'm not sure when you'd ever need this, but it's certainly an
% option. 
% 
%% Example 1: Open a figure and start clicking around to get RGB values:    
% 
% imagesc(peaks) 
% colorpicker
% 
%% Example 2: Do the same as above, but specify label format:  
% 
% imagesc(peaks) 
% colorpicker('fontweight','bold','color','r') 
% 
%% Example 3: Get 8 bit integer values:  
% 
% imagesc(peaks) 
% colorpicker(256)
% 
%% Example 4: Return an N-by-3 array of RGB values, where N is the number of times you click: 
% 
% imagesc(peaks) 
% rgb = colorpicker;
% 
%% Example 5: Get individual R, G, and B arrays scaled from 0 to 255 and specify label formatting:      
% 
% imagesc(peaks) 
% [r,g,b] = colorpicker(256,'fontangle','italic',...
%     'fontsize',30,'backgroundcolor','w'); 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% at Austin's Institute for Geophysics (UTIG), October 2015. 
% http://www.chadagreene.com. 
%
% See also ginput, ColorSpec, and impoint.  


%% Set defaults: 

rgbdivisor = 255; % Divide by 255 to match [0 1] Matlab RGB scaling. 

%% Input checks: 

if nargin>0
    if isnumeric(varargin{1})
        userrange = varargin{1}; 
        varargin(1) = []; 
        switch userrange
            case 1 
                rgbdivisor = 255; 
            case 256 
                rgbdivisor = 1; % data will be natively be 8 bit. 
            otherwise
                error('Color value range must be 1 or 256.') 
        end
    end
end

%% Get an image of the current frame: 


gca; % If no figure is open, this isinitializes a set of axes in a new figure.  

if exist('export_fig.m','file')==2
    im = export_fig(gcf,'-m1','-a1','-nocrop'); 
else
    disp('Cannot find export_fig. Using getframe instead, even though it''s a bit slower.') 
    im = getframe(gcf); 
    im = im.cdata; 
end

if ismatrix(im)==1
    R = flipud(im); 
    G = flipud(im); 
    B = flipud(im); 
else
    % Deconcatenate RGB image and flip it to match pointer units: 
    R = flipud(im(:,:,1)); 
    G = flipud(im(:,:,2)); 
    B = flipud(im(:,:,3)); 
end


%% User interface: 

% Prepare to run the UI: 
InitialNumberTitle = get(gcf,'NumberTitle'); 
InitialName = get(gcf,'Name'); 
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','ColorPicker is running. Click for RGB values or press Return to exit.') 

runui = true; 
k = 0; % initialize a counter

while runui
    w = waitforbuttonpress; 
    switch w 
        case 1 % (keyboard press) 
        	key = get(gcf,'currentcharacter'); 
            
            switch key
                case {13,27} % Return or escape keys. 
                    set(gcf,'Name',InitialName,'NumberTitle',InitialNumberTitle); 
                    try
                        delete(h)
                    end
                    runui = false; 
                    
                otherwise
                    % carry on
                    
            end
            
                
        case 0 % (mouse click)
            
            ptf = get(gcf,'CurrentPoint') ;
            pta = get(gca,'CurrentPoint') ;
            k = k+1; 
               
            try
                delete(h)
            end
            
            switch rgbdivisor
                case 255
            
                    r(k) = interp2(double(R),ptf(1),ptf(2))/255; 
                    g(k) = interp2(double(G),ptf(1),ptf(2))/255;
                    b(k) = interp2(double(B),ptf(1),ptf(2))/255;
            
                    h = text(pta(1,1),pta(1,2),['(',sprintf('%0.2f',r(k)),',',sprintf('%0.2f',g(k)),',',sprintf('%0.2f',b(k)),')'],...
                        'horiz','center','vert','bottom',varargin{:});
                
                otherwise
            
                    r(k) = interp2(R,ptf(1),ptf(2)); 
                    g(k) = interp2(G,ptf(1),ptf(2));
                    b(k) = interp2(B,ptf(1),ptf(2));
            
                    h = text(pta(1,1),pta(1,2),['(',num2str(r(k)),',',num2str(g(k)),',',num2str(b(k)),')'],...
                        'horiz','center','vert','bottom',varargin{:});
            end
    end
end

%% format output 

switch nargout
    case 0 
        clear r g b
        
    case 1 
        varargout{1} = [r(:) g(:) b(:)]; 
        
    case 3 
        varargout{1} = r(:); 
        varargout{2} = g(:); 
        varargout{3} = b(:); 
        
    otherwise
        error('Wrong number of outputs.') 
end


end

