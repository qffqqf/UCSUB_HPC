function Coordinates_orig = readcoords(filename)
fid = fopen([filename],'r+');

tline = fgetl(fid);
co_main = [];
while ischar(tline)
    aux=length(tline);
    if aux>=4
        if strcmp(tline(1:4),'GRID')                   %If line starts with Grid: new coordinate line
            if strcmp(tline(5),',')                    %if 5th character is ' ' then all coordinates on one line
                % otherwise 5th chracter is '*'
                id_last=strfind(tline,',');
                id=str2num(tline((id_last(1)+1):((id_last(2)-1))));              % take coordinate ID (is refered to in element)
                
                % Read x-coordinate
                coordx=tline((id_last(3)+1):((id_last(4)-1)));
                min_pos=strfind(coordx,'-');             % look for a '-' in the string (indicates exponential)
                min_pos=[1,min_pos];                    % add 1 in front to have at least 1 character in string ...
                
                if min_pos(end)>1&&min_pos(end)<7       % if coord-string is not filled till end
                    coordx=[coordx(1:(min_pos(end)-1)),'e',coordx((min_pos(end)):end)];  %rewrite in exponential notation and skip the last part of the string (the length of coordx is the same)
                elseif min_pos(end)>1&&min_pos(end)==7  % if coord-string is filled
                    coordx=[coordx(1:(min_pos(end)-1)),'e',coordx((min_pos(end)):end)];    % rewrite in exponential notation and delete the last siginificant digit to keep the length the same
                end
                %             Coordinates(id,1)=str2double(coordx);    % Save the x-location
                
                % Read y-coordinate
                coordy=tline((id_last(4)+1):(id_last(5)-1));
                min_pos=strfind(coordy,'-');
                min_pos=[1,min_pos];
                if min_pos(end)>1&&min_pos(end)<7
                    coordy=[coordy(1:(min_pos(end)-1)),'e',coordy((min_pos(end)):end)];
                elseif min_pos(end)>1&&min_pos(end)==7
                    coordy=[coordy(1:(min_pos(end)-1)),'e',coordy((min_pos(end)):end)];
                end
                %             Coordinates(id,2)=str2double(coordy);
                
                % Read z-coordinate
                coordz=tline((id_last(5)+1):(id_last(6)-1));
                min_pos=strfind(coordz,'-');
                min_pos=[1,min_pos];
                if min_pos(end)>1&&min_pos(end)<7
                    coordz=[coordz(1:(min_pos(end)-1)),'e',coordz((min_pos(end)):end)];         % Here not necesarry to skip last part since coordz was read to 'end' and is nog necesarrely 8 digits
                elseif min_pos(end)>1&&min_pos(end)==7
                    coordz=[coordz(1:(min_pos(end)-1)),'e',coordz((min_pos(end)):end)];
                end
                %             Coordinates(id,3)=str2double(coordz);
                
            elseif strcmp(tline(5),'*' )        % if grid-coordinates are written in 2 lines
                
                id=str2num(tline(10:17));
                coordx=tline(41:56);
                min_pos=findstr(coordx,'-');
                min_pos=[1,min_pos];
                if min_pos(end)>1
                    coordx=[coordx(1:(min_pos(end)-1)),'e',coordx((min_pos(end)):end)];
                end
                %             Coordinates(id,1)=str2double(coordx);
                
                coordy=tline(57:end);
                min_pos=findstr(coordy,'-');
                min_pos=[1,min_pos];
                if min_pos(end)>1
                    coordy=[coordy(1:(min_pos(end)-1)),'e',coordy((min_pos(end)):end)];
                end
                %             Coordinates(id,2)=str2double(coordy);
                
                tline = fgetl(fid);
                
                coordz=tline(9:end);
                min_pos=findstr(coordz,'-');
                min_pos=[1,min_pos];
                if min_pos(end)>1
                    coordz=[coordz(1:(min_pos(end)-1)),'e',coordz((min_pos(end)):end)];
                end
                %             Coordinates(id,3)=str2double(coordz);
            end
            
            co_main = [co_main; [id, str2double(coordx), str2double(coordy), str2double(coordz)]];
            
        end
    end
    tline = fgetl(fid);
end

fclose(fid);

Coordinates_orig= co_main;

% clear fid
% clear id
% clear tline
% clear ans
% clear Coordinates

% save Coords_orig_leg3f_materialise_9_29.mat