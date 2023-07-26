# How to use myQAMmod.c and myQAMdemod.c

# Use the myQAMmod MEX functions

>> modulatedSignal = myQAMmod([0:3], 4); % This will call your qammod MEX function
>> matlabSignal = qammod([0:3], 4, 'gray'); % This will call MATLAB qammod function
>> difference = sum(modulatedSignal - matlabSignal) % Result should be 0

# Use the myQAMdemod MEX functions

>> data = [0:3];
>> modData = qammod(data, 4, 'gray'); % This will call MATLAB qammod function
>> demodData = myQAMdemod(modData, 4); % This will call your qamdemod MEX function to demod modData
>> difference = sum(data) - sum(demodData) % Result should be 0
