import matplotlib.pyplot as plt

# Read the data from the files
with open('mach.dat', 'r') as mach_file:
    mach_data = [float(line.strip()) for line in mach_file]

with open('pressure.dat', 'r') as pressure_file:
    pressure_data = [float(line.strip()) for line in pressure_file]

with open('cell_centers.dat', 'r') as cell_centers_file:
    cell_centers_data = [float(line.strip()) for line in cell_centers_file]

# Plot the data
plt.plot(cell_centers_data, mach_data, label='Mach')
plt.plot(cell_centers_data, pressure_data, label='Pressure')

# Add labels and legend
plt.xlabel('Cell Centers')
plt.ylabel('Values')
plt.legend()

# Show the plot
plt.show()
