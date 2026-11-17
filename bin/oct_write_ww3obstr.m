function [messg,errno] = oct_write_ww3obstr(fname,d1,d2)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  23-Oct-2012 |                      |
%|                    +----------------------------+                      | 
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| Write the output arrays into ascii file                                |
%|                                                                        |
%| oct_write_ww3file(fname,d1,d2)                                             |
%|                                                                        |
%| INPUT                                                                  |
%|  fname       : Output file name                                        |
%|  d1,d2       : Output 2D obstruction arrays in x (d1) and y (d2)       |
%|                                                                        |
%| OUTPUT                                                                 |
%|  messg       : Error message. Is blank if no error occurs              |
%|  errno       : Error number. Is zero for succesful write               |
% -------------------------------------------------------------------------

[Ny,Nx] = size(d1);

fid = fopen(fname,'w');

[messg,errno] = ferror(fid);

if (errno == 0)
   %@@@ one fprintf call per array (same text as row by row)
   fmt = [repmat(' %d ',1,Nx) '\n'];
   fprintf(fid,fmt,d1.');
   fprintf(fid,'\n');
   fprintf(fid,fmt,d2.');
else
   fprintf(1,'!!ERROR!!: %s \n',messg);
end;

fclose(fid);

return;
