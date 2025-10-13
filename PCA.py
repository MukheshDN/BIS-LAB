import numpy as np

def parallel_cellular_routing(grid_size=10, alpha=0.5, max_iter=100):
    rows, cols = grid_size, grid_size
    cost = np.full((rows, cols), 9999.0)   # finite large initial cost
    new_cost = np.full((rows, cols), 9999.0)
    link_cost = 1                          # cost to move between neighbors

    source = (0, 0)
    destination = (rows - 1, cols - 1)
    cost[destination] = 0                  # destination cost = 0

    neighbors = [(-1,0), (1,0), (0,-1), (0,1)]

    # Parallel update loop
    for iteration in range(max_iter):
        for i in range(rows):
            for j in range(cols):
                if (i, j) == destination:
                    new_cost[i, j] = 0
                    continue

                neighbor_costs = []
                for di, dj in neighbors:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols:
                        neighbor_costs.append(cost[ni, nj] + link_cost)

                if neighbor_costs:
                    best_neighbor_cost = min(neighbor_costs)
                    new_cost[i, j] = cost[i, j] + alpha * (best_neighbor_cost - cost[i, j])

        cost = new_cost.copy()

    # Path reconstruction
    path = [source]
    current = source
    while current != destination:
        i, j = current
        next_cell = None
        min_val = float('inf')

        for di, dj in neighbors:
            ni, nj = i + di, j + dj
            if 0 <= ni < rows and 0 <= nj < cols and cost[ni, nj] < min_val:
                min_val = cost[ni, nj]
                next_cell = (ni, nj)

        if next_cell is None or next_cell == current:
            break  # safety break
        path.append(next_cell)
        current = next_cell

    return path


# Run the algorithm
best_path = parallel_cellular_routing(grid_size=10, alpha=0.5, max_iter=100)
print("Optimal Path (from source to destination):")
print(best_path)
