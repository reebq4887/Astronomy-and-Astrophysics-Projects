import glob

path = 'Processed Kepler Data/lightcurves/short cadence'
files = glob.glob(f'{path}/* Q8.txt', recursive=True)



lines = []
print(len(lines))
for file in files:
    with open(file, 'r') as f:
        lines.append(f.readlines()[1:])        

for line in lines:
    with open(f'{path}/combined_Q8.txt', 'a') as f:
        f.writelines(line)