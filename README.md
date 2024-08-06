# ACOvsBruteForce
In this project we are trying to find the shortest path by using 2 different methods (ACO and Brute Force Methods)

Here's the README file explanation in English:

---

# Migros Delivery Optimization

## Description

In this project, we aim to find the shortest path in the **Migros Delivery** application using two different methods: **Ant Colony Optimization (ACO)** and **Brute Force**. The project utilizes these algorithms to determine the shortest route between various coordinate points.

## Features

- **Coordinate Reading**: Coordinates are read from a file and imported into the application.
- **Shortest Path Calculation**: The user can select between two different optimization methods to find the shortest route:
  - **Brute Force**: Tests all possible paths to find the shortest one.
  - **Ant Colony Optimization**: Utilizes pheromones to explore possible paths.
- **Visualization**: The found routes are drawn on a graph, allowing users to visually inspect the results.

## Usage

1. Specify the coordinates in a file named `input02.txt`. Each line should contain a coordinate pair (x,y).
2. Run the application:
   ```
   java BaranOnder
   ```
3. The user can select the optimization method to find the shortest path (the default method is set to 2 in the code).

## Algorithms

### Brute Force Method
This method calculates all possible routes to find the shortest distance. When a shorter distance is found, the application presents the shortest path and distance information to the user.

### Ant Colony Optimization (ACO)
This method is based on the principle of virtual ants leaving pheromones and following these pheromones to find the best path. In each iteration, ants discover new paths, and pheromone levels are updated.
