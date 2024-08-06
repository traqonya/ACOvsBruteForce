import java.awt.*;
import java.util.*;
import java.io.*;
import java.nio.file.*;
import java.util.List;
/**
 * @author baranonder , Student ID: 2022400111
 * @version 1.0
 * @since Date: 13.05.2024
 * Migros Delivery: In this assignment we are trying to
 * find the shortest path by using 2 different methods (ACO and
 * Brute Force Methods).
*/
class Coordinate {
    double xPos, yPos;

    public Coordinate(double x, double y) {
        this.xPos = x; // Set x coordinate
        this.yPos = y; // Set y coordinate
    }

    // Returns the distance to another point
    public double measureDistance(Coordinate other) {
        return Math.sqrt(Math.pow(this.xPos - other.xPos, 2) + Math.pow(this.yPos - other.yPos, 2));
    }

    // Computes the end point of a line from an angle and radius
    public Coordinate calculateBorder(double angle, double radius) {
        return new Coordinate(xPos + radius * Math.cos(angle), yPos + radius * Math.sin(angle));
    }
}

public class BaranOnder {
    static final int MAX_ITERATIONS = 100; // Max number of optimization cycles
    static final int NUM_ANTS = 50; // Number of ants in the simulation
    static final double EVAPORATION_RATE = 0.9; // Rate of pheromone evaporation
    static final double ALPHA = 0.8; // Influence of pheromone in path selection
    static final double BETA = 5.5; // Influence of heuristic information
    static final double INITIAL_PHEROMONE = 0.1; // Initial pheromone level on paths
    static final double Q = 0.0001; // Pheromone deposit quantity
    static final double RADIUS = 0.01; // Size of the dots representing the coordinates
    static int chosenMethod = 2; // Selected optimization method
    static int antColonyGraph = 2; // Graph mode for visualization
    private List<Coordinate> points = new ArrayList<>(); // List of all points
    private double[][] pheromones; // Pheromone levels between points
    private List<Integer> optimalPath; // Best path found
    private double optimalDistance = Double.MAX_VALUE; // Shortest distance found
    private Random random = new Random(); // Random generator for probability calculations

    public BaranOnder() {
        this.pheromones = new double[30][30]; // Initialize pheromone matrix
        StdDraw.setCanvasSize(800, 800); // Set the size of the drawing canvas
        StdDraw.setXscale(0, 1); // Set the x scale of the canvas
        StdDraw.setYscale(0, 1); // Set the y scale of the canvas
    }

    // Reads coordinates from a file and initializes points list
    public void loadCoordinates(String filePath) throws IOException {
        List<String> lines = Files.readAllLines(Paths.get(filePath)); // Read all lines from the file
        for (String line : lines) {
            String[] coordinates = line.split(","); // Split line into x and y parts
            double x = Double.parseDouble(coordinates[0].trim()); // Parse x coordinate
            double y = Double.parseDouble(coordinates[1].trim()); // Parse y coordinate
            points.add(new Coordinate(x, y)); // Add new point to the list
        }
        setupPheromones(); // Initialize pheromones after loading points
    }

    // Fills the pheromone matrix with initial values
    private void setupPheromones() {
        for (int i = 0; i < points.size(); i++) {
            Arrays.fill(pheromones[i], INITIAL_PHEROMONE); // Set initial pheromone level for each path
        }
    }

    // Selects and executes the optimization routine based on the chosen method
    public void optimizeRoute(int method) {
        if (method == 1) {
            calculateShortestRoute(); // Brute-force approach
        } else if (method == 2) {
            simulateAntColony(); // Ant colony optimization
        }
    }

    // Brute-force algorithm to find the shortest possible route
    private void calculateShortestRoute() {
        long start = System.currentTimeMillis(); // Record start time
        int n = points.size(); // Number of points
        int[] indices = new int[n - 1]; // Array to hold indices for permutation
        for (int i = 1; i < n; i++) {
            indices[i - 1] = i; // Initialize indices array
        }

        optimalDistance = Double.MAX_VALUE; // Reset the shortest distance
        do {
            double currentDistance = calculateTotalDistance(indices); // Calculate distance for current permutation
            if (currentDistance < optimalDistance) {
                optimalDistance = currentDistance; // Update shortest distance
                optimalPath = new ArrayList<>(); // Create new list for best path
                optimalPath.add(0); // Add starting point
                for (int index : indices) optimalPath.add(index); // Add points in current order
                optimalPath.add(0); // Complete the loop to the start
            }
        } while (nextPermutation(indices)); // Generate next permutation

        long end = System.currentTimeMillis(); // Record end time
        System.out.println("Method: Brute-Force");
        System.out.printf("Shortest Distance: %.5f\n", optimalDistance);
        System.out.printf("Shortest Path: %s\n", formatRoute(optimalPath));
        System.out.println("Time it takes to find the shortest path: " + ((end - start) / 1000.0) + " seconds");

        drawRoute(optimalPath, true); // Draw the best route found
    }

    // Formats the route for display, converting indices to 1-based
    private String formatRoute(List<Integer> route) {
        List<Integer> formatted = new ArrayList<>();
        for (int i : route) formatted.add(i + 1); // Convert to 1-based index
        return formatted.toString();
    }

    // Generates the next permutation of indices in-place
    private boolean nextPermutation(int[] array) {
        int i = array.length - 2; // Start from the second last element
        while (i >= 0 && array[i] >= array[i + 1]) i--; // Find the first decreasing element
        if (i == -1) return false; // No more permutations

        int j = array.length - 1; // Start from the last element
        while (array[j] <= array[i]) j--; // Find the element to swap with
        swap(array, i, j); // Swap them
        reverse(array, i + 1, array.length - 1); // Reverse the sequence after i
        return true; // Permutation generated
    }

    // Swaps two elements in an array
    private void swap(int[] array, int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    // Reverses a portion of an array
    private void reverse(int[] array, int start, int end) {
        while (start < end) swap(array, start++, end--);
    }

    // Calculates the total distance of a route starting at the first index
    private double calculateTotalDistance(int[] indices) {
        double length = points.get(0).measureDistance(points.get(indices[0])); // Start with the first leg of the journey
        for (int i = 0; i < indices.length - 1; i++) {
            length += points.get(indices[i]).measureDistance(points.get(indices[i + 1])); // Add distance for each leg
            if (length >= optimalDistance) return Double.MAX_VALUE; // If the current length exceeds the best, abort
        }
        length += points.get(indices[indices.length - 1]).measureDistance(points.get(0)); // Complete the loop
        return length; // Return the total length
    }

    // Simulates ant colony optimization to find the shortest route
    private void simulateAntColony() {
        long start = System.currentTimeMillis(); // Start timing
        for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) { // Loop for a fixed number of iterations
            for (int ant = 0; ant < NUM_ANTS; ant++) { // Simulate each ant
                List<Integer> route = buildRouteForAnt(); // Build a route for this ant
                double routeLength = calculateTotalDistance(toIntArray(route)); // Calculate the length of the route
                if (routeLength < optimalDistance) {
                    optimalDistance = routeLength; // Update the best distance
                    optimalPath = new ArrayList<>(route); // Update the best path
                }
                updatePheromones(route, routeLength); // Update pheromones based on this route
            }
            evaporatePheromones(); // Reduce all pheromone levels
        }
        long end = System.currentTimeMillis(); // End timing
        double duration = (end - start) / 1000.0; // Calculate duration
        System.out.println("Method: Ant-Colony Method");
        System.out.println("Shortest Distance: " + optimalDistance);
        System.out.println("Shortest Path: " + formatRoute(optimalPath));
        System.out.println("Time it takes to find the shortest path: " + duration + " seconds");

        // Draw pheromone levels if required
        if (antColonyGraph == 2) {
            drawPheromoneLevels(); // Draw only pheromone levels
        }
    }

    // Converts a list of integers to an array of integers
    private int[] toIntArray(List<Integer> list) {
        return list.stream().mapToInt(i -> i).toArray(); // Convert List<Integer> to int[]
    }

    // Builds a route for an ant based on pheromone levels and heuristic information
    private List<Integer> buildRouteForAnt() {
        List<Integer> route = new ArrayList<>(); // Initialize the route list
        boolean[] visited = new boolean[points.size()]; // Track visited points
        int current = 0; // Start at the first point
        route.add(current); // Add start point to route
        visited[current] = true; // Mark start point as visited

        while (route.size() < points.size()) { // Until all points are visited
            int next = selectNextPoint(current, visited); // Select the next city based on pheromones and distance
            if (next == -1) break; // If no valid next city, break loop
            route.add(next); // Add next city to route
            visited[next] = true; // Mark this city as visited
            current = next; // Move to the next city
        }
        route.add(0); // Return to the start point at the end
        return route; // Return the completed route
    }

    // Selects the next point for an ant to visit
    private int selectNextPoint(int current, boolean[] visited) {
        double total = 0.0; // Total sum of probabilities
        double[] probabilities = new double[points.size()]; // Probabilities of moving to each point
        for (int j = 0; j < points.size(); j++) {
            if (!visited[j]) {
                probabilities[j] = Math.pow(pheromones[current][j], ALPHA) * Math.pow(1.0 / points.get(current).measureDistance(points.get(j)), BETA); // Calculate the probability based on pheromone level and heuristic
                total += probabilities[j]; // Add to total probability
            } else {
                probabilities[j] = 0; // If already visited, probability is zero
            }
        }

        double r = random.nextDouble() * total; // Randomly select a probability threshold
        double sum = 0.0; // Cumulative probability
        for (int j = 0; j < probabilities.length; j++) {
            sum += probabilities[j]; // Add this point's probability
            if (sum > r) {
                return j; // If cumulative probability exceeds threshold, select this point
            }
        }
        return -1; // If no point selected, return -1 (should not happen)
    }

    // Updates pheromone levels based on a given route and its length
    private void updatePheromones(List<Integer> route, double routeLength) {
        double deposit = Q / routeLength; // Calculate pheromone deposit
        for (int i = 0; i < route.size() - 1; i++) {
            int u = route.get(i);
            int v = route.get(i + 1);
            pheromones[u][v] += deposit; // Deposit pheromones on this path
            pheromones[v][u] += deposit; // Deposit pheromones on reverse path
        }
    }

    // Reduces pheromone levels on all paths
    private void evaporatePheromones() {
        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points.size(); j++) {
                pheromones[i][j] *= EVAPORATION_RATE; // Evaporate pheromones
            }
        }
    }

    // Draws the route based on the specified method
    private void drawRoute(List<Integer> route, boolean isBruteForce) {
        if (isBruteForce) {
            drawRouteBruteForce(route); // Draw using brute force method
        } else {
            drawRouteAntColony(route); // Draw using ant colony method
        }
        StdDraw.show(); // Display the canvas
    }

    // Draws the route found by brute force
    private void drawRouteBruteForce(List<Integer> route) {
        StdDraw.clear(); // Clear the canvas
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        drawRoutes(route); // Draw all lines between points

        // Draw all points
        for (int i = 0; i < route.size(); i++) {
            int current = route.get(i);
            Coordinate currentPoint = points.get(current);
            Color pointColor = StdDraw.GRAY; // Default color for all points
            drawPoint(currentPoint, String.valueOf(current + 1), pointColor); // Draw each point
        }

        // Redraw the starting point to ensure it is highlighted
        Coordinate startPoint = points.get(route.get(0)); // Get the starting point
        drawPoint(startPoint, "1", StdDraw.PRINCETON_ORANGE); // Highlight the starting point

        StdDraw.show(); // Display the canvas
    }

    // Draws the route found by ant colony optimization
    private void drawRouteAntColony(List<Integer> route) {
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        for (int i = 1; i < route.size(); i++) {
            Coordinate current = points.get(route.get(i)); // Get current point
            Coordinate previous = points.get(route.get(i - 1)); // Get previous point
            drawLine(previous, current, RADIUS, 0.01, StdDraw.BLACK); // Draw line between them
        }

        // Ensure points are drawn after lines
        for (Coordinate point : points) {
            drawPoint(point, String.valueOf(points.indexOf(point) + 1), StdDraw.GRAY); // Draw each point
        }
        StdDraw.show(); // Display the canvas
    }

    // Draws all lines between points in a route
    private void drawRoutes(List<Integer> route) {
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        StdDraw.setPenColor(StdDraw.BLACK); // Set line color to black
        StdDraw.setPenRadius(0.005); // Set line thickness

        for (int i = 1; i < route.size(); i++) {
            Coordinate currentPoint = points.get(route.get(i)); // Get current point
            Coordinate previousPoint = points.get(route.get(i - 1)); // Get previous point
            drawLine(previousPoint, currentPoint, RADIUS, 0.005, StdDraw.BLACK); // Draw line between them
        }

        // Optionally complete the route loop, but do NOT draw points here
        if (route.size() > 1) {
            Coordinate firstPoint = points.get(route.get(0)); // Get the first point
            Coordinate lastPoint = points.get(route.get(route.size() - 1)); // Get the last point
            drawLine(lastPoint, firstPoint, RADIUS, 0.005, StdDraw.BLACK); // Draw line to complete the loop
        }
    }

    // Finds the maximum pheromone level on any path
    private double findMaxPheromone() {
        double max = 0; // Initialize max to zero
        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points.size(); j++) {
                if (pheromones[i][j] > max) {
                    max = pheromones[i][j]; // Update max if a higher level is found
                }
            }
        }
        return max; // Return the maximum level found
    }

    // Draws pheromone levels between points
    private void drawPheromoneLevels() {
        StdDraw.clear(); // Clear the canvas
        double maxPheromone = findMaxPheromone(); // Find the maximum pheromone level
        if (maxPheromone == 0) return; // If no pheromones, exit

        StdDraw.enableDoubleBuffering(); // Enable double buffering
        for (int i = 0; i < points.size(); i++) {
            for (int j = i + 1; j < points.size(); j++) {
                double pheromoneLevel = pheromones[i][j]; // Get pheromone level between two points
                if (pheromoneLevel > 0) {
                    Coordinate from = points.get(i); // Get starting point
                    Coordinate to = points.get(j); // Get ending point
                    double normalizedPheromone = pheromoneLevel / maxPheromone; // Normalize pheromone level
                    double lineThickness = 0.0001 + 0.002 * normalizedPheromone; // Calculate line thickness
                    int alpha = (int) (255 * normalizedPheromone); // Calculate transparency
                    Color pheromoneColor = new Color(0, 0, 0, Math.min(255, alpha)); // Create color with transparency
                    drawLine(from, to, lineThickness, pheromoneColor); // Draw pheromone line
                }
            }
        }

        // Draw points on top of pheromone levels
        drawPoints(); // Draw all points

        StdDraw.show(); // Display the canvas
    }

    // Draws all points on the canvas
    private void drawPoints() {
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        for (Coordinate point : points) {
            drawPoint(point, String.valueOf(points.indexOf(point) + 1), StdDraw.GRAY); // Draw each point with its index
        }
        StdDraw.show(); // Display the canvas
        // Optionally highlight the starting point (assuming the first point is the start)
    }

    // Draws a line between two points with specified thickness and color
    private void drawLine(Coordinate from, Coordinate to, double thickness, Color color) {
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        StdDraw.setPenColor(color); // Set pen color
        StdDraw.setPenRadius(thickness); // Set pen thickness
        StdDraw.line(from.xPos, from.yPos, to.xPos, to.yPos); // Draw the line
    }

    // Draws a point with a label and specified color
    private void drawPoint(Coordinate point, String label, Color color) {
        StdDraw.setPenColor(color); // Set pen color for the point
        StdDraw.filledCircle(point.xPos, point.yPos, RADIUS); // Draw filled circle for the point
        StdDraw.setPenColor(StdDraw.BLACK); // Set pen color for the label
        StdDraw.text(point.xPos, point.yPos, label); // Draw the label at the point
    }

    // Draws a line with specified radius and thickness using an angle calculation
    private void drawLine(Coordinate from, Coordinate to, double radius, double thickness, Color color) {
        StdDraw.enableDoubleBuffering(); // Enable double buffering
        StdDraw.setPenRadius(thickness); // Set pen thickness
        StdDraw.setPenColor(color); // Set pen color
        double angle = Math.atan2(to.yPos - from.yPos, to.xPos - from.xPos); // Calculate the angle
        Coordinate start = from.calculateBorder(angle, radius); // Calculate start point
        Coordinate end = to.calculateBorder(angle + Math.PI, radius); // Calculate end point
        StdDraw.line(start.xPos, start.yPos, end.xPos, end.yPos); // Draw the line
    }

    public static void main(String[] args) {
        try {
            BaranOnder optimizer = new BaranOnder();
            optimizer.loadCoordinates("input02.txt"); // Load coordinates from file

            // Execute the optimization routine based on selected method
            optimizer.optimizeRoute(chosenMethod);

            // Conditionally draw the route based on the graph setting
            if (chosenMethod == 2 && antColonyGraph != 2) {
                optimizer.drawRoute(optimizer.optimalPath, false); // Draw the route if not displaying pheromone levels
            }
        } catch (IOException e) {
            System.out.println("Error reading the coordinates file: " + e.getMessage()); // Handle potential IO errors
        }
    }
}
