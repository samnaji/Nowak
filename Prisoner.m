Columns = [1 2
    3 4];
Coloring_map = [ 0 0 1   
        0 1 0
        1 1 0
        1 0 0];
 
End_time = 100;b_parameter= 1.77;SIZE = 100;TotalNumberOfCells = SIZE*SIZE;
p = 0.10; Payoff = zeros(2,2); Payoff(1,1) = 1; Payoff(2,1) = b_parameter;
Transcation = zeros(SIZE,SIZE);
NE = zeros(SIZE,SIZE);
 
move = zeros(SIZE,SIZE);
move_cube = zeros(SIZE,SIZE,End_time);
 
 
% For initial conditions has two opions 
E = ones(SIZE,SIZE);A = rand(SIZE,SIZE);
 
%1-Randomly located
% % Changing E according to the initial proportion of defectors
 I = find(A < p);
 E(I) = 2;
 
%2-Fixed at the center
%E(SIZE/2,SIZE/2) = 2;

%3-Fixed at the corners

% E(1,1) = 2;
% E(1,SIZE) = 2;
% 
% E(SIZE,1) = 2;
% E(SIZE,SIZE) = 2;
for genneration_i=1:End_time
    fprintf('Step:%i\n ',[genneration_i]);
 
    for i=1:SIZE
        for j=1:SIZE
            pa = 0;
            for k=-1:1
                for h=-1:1
                    a = findingBoundary(i+k,SIZE); b = findingBoundary(j+h,SIZE);
                    pa = pa + Payoff(E(i,j), E(a,b));
                end
            end
            Transcation(i,j) = pa;
        end
    end
    for i=1:SIZE
        for j=1:SIZE
            pay = Transcation(i,j);
            NE(i,j) = E(i,j);
            for k=-1:1
                for h=-1:1
                    a = findingBoundary(i+k,SIZE); b = findingBoundary(j+h,SIZE);
                    % If the neighbour performed better
                    if (Transcation(a,b) > pay)
                        pay = Transcation(a,b);
                        NE(i,j) = E(a,b);
                    end
                end
            end
        end
    end
    for i=1:SIZE
        for j=1:SIZE
            move(i,j) = Columns(NE(i,j),E(i,j));
            E(i,j) = NE(i,j);
        end
    end
    move_cube(:,:,genneration_i) = move;
    
 
end 
 %% gif
for genneration_i=1:End_time
    imagesc(move_cube(:,:,genneration_i));
    colormap(Coloring_map)
    title(['Generation ',num2str(genneration_i)])
    axis('off')
    set(gcf,'color','w');
%     if  max(genneration_i==[1  20 50 100])
%     saveas(gca,['Fig3A_Generation_' num2str(genneration_i) '.png'])
%     end
      F(genneration_i) = getframe;
    % Capture the plot as an image 
      axis('off')
      set(gcf,'color','w');
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      filename= [num2str(SIZE) '_b_' num2str(100*b_parameter) '_p' num2str(100*p)  '_rand.gif'];
    
      if genneration_i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
 
end
 
function y = findingBoundary(cell,size)
 
switch cell
    case 0 
        y = size;
    case size+1 
        y = 1; 
    otherwise
        y = cell;
end
end
