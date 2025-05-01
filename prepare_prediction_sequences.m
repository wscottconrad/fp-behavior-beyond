% this function prepares sequences/windows of data
% For predictMovementfromFluorescence
% W Scott Conrad 14/2/25 <3 made with Claude


function [X, Y] = prepare_prediction_sequences(fluorescence, speed, ...
    session_indices, sequence_length, prediction_horizon)
    % X will contain sequences of fluorescence data
    % Y will contain future movement data we want to predict
    
    % Initialize arrays
    X = cell(1, length(session_indices));  % Input sequences
    Y = cell(1, length(session_indices));  % Output sequences (speed and angular velocity)
    
    for i = 1:length(session_indices)
        
        fluor_session = fluorescence(i,:);
        if rem(i,2) ~=0
            speed_session = speed(i/2+0.5,:);
        else
            speed_session = speed(i/2,:);
        end
%         angular_vel_session = angular_vel(session_mask);
        
        
        % Calculate how many complete sequences we can make
        num_sequences = length(fluor_session) - sequence_length - sequence_length + 1;
        
        % Initialize arrays
        X{1, i} = zeros(sequence_length, num_sequences);  % Input sequences
        Y{1, i} = zeros(sequence_length, num_sequences);  % Output sequences (speed)
%         Y{end+1} = zeros(2, prediction_horizon, num_sequences);  % Output sequences (speed and angular velocity)

        
        if num_sequences > 0
            for j = 1:num_sequences
                % Create sequence from this session
                X{1,i}(:, j) = fluor_session(j:j+sequence_length-1);
                Y{1,i}(:, j) = speed_session(j+sequence_length:j+sequence_length+sequence_length-1);
                    % Y{end}(2, :, j) = angular_vel_session(j+sequence_length:j+sequence_length+prediction_horizon-1);
               
                
            end
        end
    end
    
    % Convert cell arrays to matrices for training
    X = cat(2, X{:});
    Y = cat(2, Y{:}); 
    
%     % Create sequences
%     for i = 1:num_sequences
%         % Input: window of fluorescence data
%         X(1, :, i) = fluorescence(i:i+sequence_length-1);
%         
%         % Output: future movement data
%         movement_start = i + sequence_length;
%         movement_end = movement_start + prediction_horizon - 1;
%         Y(1, :, i) = speed(movement_start:movement_end);
%         Y(2, :, i) = angular_vel(movement_start:movement_end);
%     end
end