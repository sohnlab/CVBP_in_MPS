function result = simulateModel(spicePath, fileName, filePath, varargin)
% SIMULATEMODEL is used to simulate a spice netlist. It takes the following
% Inputs: 
% spicePath: The path to the LT SPice installation start file
% for example 'C:\Program Files\LTC\LTspiceXVII\start.exe' This should be a
% string
% filename: Name of the netlist file to be simulated. Do not use .net
% extension. This can be a string / char or even an int or float
% filepath: Provide the complete file path to the netlist and end it with a
% backslash (\). This has to be a string.

if ischar(fileName)
    filename = fileName;
else
    filename = num2str(fileName);
end

outputFile = sprintf('%s%s.raw', filePath, filename);

if ispc
    command = sprintf('"%s" -b "%s%s.net"', spicePath, filePath, filename);
elseif isunix
    command = sprintf('cd "%s" && wine "%s" -Run -b "%s.net"', filePath, spicePath, filename);
end
system(command, '-echo'); % run the command

while ~exist(outputFile, 'file')
    pause(0.01);
end
result = LTspice2Matlab(outputFile, varargin{1});

% clean up files
if ispc
    command = sprintf('del "%s%s.*"', filePath, filename);
elseif isunix
    command = sprintf('cd "%s" && rm %s.*', filePath, filename);
end
system(command, '-echo');

end
