import numpy as np

# Input the endpoints
L1 = 0
L2 = 7

# Number of divisions (elements)
num_ele = 10

# Linear(1) or Quad(2)
choice = 1

# [Takes effect for Quad only] Determine the location of the midpoint, between [-0.99, 0.99]
midpoint_loc = 0.1 # A little more to the right

# Produce the array
if choice == 1: # Linear Generation
    output_mesh = np.linspace(L1, L2, num_ele + 1)
elif choice == 2: # Quadratic Generation
    inter_mesh = np.linspace(L1, L2, num_ele)
    output_mesh = [L1]
    for i in range(1, num_ele):
        p1, p2 = inter_mesh[i], inter_mesh[i - 1]
        mid = (p1 + p2)/2
        # Perform the sliding
        mid *= 1 + midpoint_loc
        output_mesh.extend([mid, p1])

    output_mesh = np.array(output_mesh)

# Write in basic mesh format
filename = "input_mesh.in"
with open(filename, "w") as m:
    # Write the number of nodes
    m.write(f"{len(output_mesh)}\n")
    for node in output_mesh:
        m.write(f"{node}\n")

