function [] = plotTracks_Color(Tracks,Color,LineWidth);

%%Input Tracks as a cell with a list of trajectories that must be Nx3,
%%with columns being: t, x, y


n_tracks = size(Tracks,1);

    % 2D case
    for i = 1 : n_tracks;
        
        track = Tracks{i};
        
        x = track(:,2);
        y = track(:,3);
        
        plot(x, y, 'Color', Color,'LineWidth',LineWidth);
        hold on
        
    end
    


end