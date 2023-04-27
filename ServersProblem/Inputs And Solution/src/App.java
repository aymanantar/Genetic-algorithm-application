import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class App {
    static int V, E, R, C, X;
    static int[] vidSize;
    static int[] dataCenterCost;
    static int[][] serverCost;
    static int[][] requestFrequency;
    static int[][] CellsScore;
    public static void main(String[] args) {
        Read();
        CellsScore = CalculateCellScores();
        // Hill Climbing method
        // boolean [][] Sol = hillClimbing();
        // Genatic algorithm method
        boolean [][] Sol= getLeastScoreIndividual(GeneticAlgorithm());
        System.out.println( Scoringfunction(Sol) +" is the highest score achieved, the answer is as follows: ");
        for (int j = 0 ; j<C; j++)
        {
            for (int k =0 ; k< V ; k++)
            {
                if (Sol[j][k]) {
                  //  System.out.print("1 ");
                    }
                    else {
                     //   System.out.print("0 ");
                    }
            }
           // System.out.println("  ");
        }
    }
    public static int[][] CalculateCellScores() {
        int numEndpoints = requestFrequency.length;
        int numVideos = requestFrequency[0].length;
        int numCaches = serverCost[0].length; // Assumes serverCost is a square matrix
        
        int[][] cellScores = new int[numCaches][numVideos]; // Matrix to store cell scores
        
        int[] bestLatency = new int[numEndpoints]; // Store best latencies for each endpoint
        
        // Pre-process cache latencies for each endpoint
        for (int endpointId = 0; endpointId < numEndpoints; endpointId++) {
            bestLatency[endpointId] = dataCenterCost[endpointId]; // Initialize with data center cost
            for (int cacheId = 0; cacheId < numCaches; cacheId++) {
                for (int vidId = 0; vidId < numVideos; vidId++) {
                    if (serverCost[endpointId][cacheId] < bestLatency[endpointId]) {
                        bestLatency[endpointId] = serverCost[endpointId][cacheId];
                    }
                }
            }
        }
    
        // Calculate cell scores
        for (int endpointId = 0; endpointId < numEndpoints; endpointId++) {
            for (int vidId = 0; vidId < numVideos; vidId++) {
                int currentRequestFrequency = requestFrequency[endpointId][vidId];
                double currentScore = (currentRequestFrequency * (dataCenterCost[endpointId] - bestLatency[endpointId]));
                for (int cacheId = 0; cacheId < numCaches; cacheId++) {
                    // Update cell score with estimated effect on solution score
                    cellScores[cacheId][vidId] += currentScore;
                }
            }
        }
    
        return cellScores;
    }    
    public static boolean[][] getLeastScoreIndividual(boolean[][][] population) {
        double maxScore = Double.MIN_VALUE;
        boolean[][] leastScoreIndividual = null;
    
        for (boolean[][] individual : population) {
            // Evaluate the fitness score for each individual
            double score = Scoringfunction(individual);
    
            // Update the minimum score and least score individual if necessary
            if (score > maxScore) {
                maxScore = score;
                leastScoreIndividual = individual;
            }
        }
    
        return leastScoreIndividual;
    }
    public static void Read() {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        System.out.println("Enter the input file name:");
        String fileName = "";
        try {
            fileName = "inputs/"+ br.readLine() + ".txt";
        } catch (IOException e) {
            e.printStackTrace();
        }
        try (BufferedReader fileReader = new BufferedReader(new FileReader(fileName))) {
            String line = fileReader.readLine();
            String[] parts = line.split(" ");
            V = Integer.parseInt(parts[0]);
            E = Integer.parseInt(parts[1]);
            R = Integer.parseInt(parts[2]);
            C = Integer.parseInt(parts[3]);
            X = Integer.parseInt(parts[4]);

            vidSize = new int[V];
            dataCenterCost = new int[E];
            serverCost = new int[E][C];
            requestFrequency = new int[E][V];
            for (int i=0 ; i<E; i++)
            {
                for (int j=0; j<C ; j++)
                {
                    serverCost[i][j]=100000000;
                }
            }
            for (int i=0 ; i<E; i++)
            {
                for (int j=0; j<V ; j++)
                {
                    requestFrequency[i][j]=0;
                }
            }
            line = fileReader.readLine();
            parts = line.split(" ");
            for (int i = 0; i < V; i++) {
                vidSize[i] = Integer.parseInt(parts[i]);
            }

            for (int i = 0; i < E; i++) {
                line = fileReader.readLine();
                parts = line.split(" ");
                dataCenterCost[i] = Integer.parseInt(parts[0]);
                int connectedServers = Integer.parseInt(parts[1]);
                for (int j = 0; j < connectedServers; j++) {
                    line = fileReader.readLine();
                    parts = line.split(" ");
                    int id = Integer.parseInt(parts[0]);
                    int cost = Integer.parseInt(parts[1]);
                    serverCost[i][id] = cost;
                }
            }

            for (int i = 0; i < R; i++) {
                line = fileReader.readLine();
                parts = line.split(" ");
                int vidId = Integer.parseInt(parts[0]);
                int endPointId = Integer.parseInt(parts[1]);
                int requests = Integer.parseInt(parts[2]);
                requestFrequency[endPointId][vidId] = requests;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public static boolean isValid(boolean[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
    
        for (int i = 0; i < rows; i++) {
            int Count = 0;
            for (int j = 0; j < cols; j++) {
                if (matrix[i][j]) {
                    Count+=vidSize[j];
                }
            }
            if (Count > X) {
                return false; // If true count in a row exceeds maxTrueValues, return false
            }
        }
    
        return true; // If all rows pass the condition, return true
    }
    public static boolean[][][] GeneticAlgorithm() {
        int populationSize = 10;
        int maxGenerations = 10;
        double crossoverProb = 5;
        double mutationProb = 50;
        double minFitnessScoreRatioThreshold =0.001;
        double oldmax = 1; 
        boolean[][][] population = new boolean[populationSize][C][V];
        Random random = new Random();
        // Initialize the population randomly
        for (int i = 0; i < populationSize; i++) {
            for (int j = 0; j < C; j++) {
                for (int k = 0; k < V; k++) {
                    population[i][j][k] = false;
                }
            }
        }
        for (int generation = 1; generation <= maxGenerations; generation++) {
            List<boolean[][]> newPopulation = new ArrayList<>();
            for (int  i= 0 ; i<populationSize ; i++){newPopulation.add(population[i]);}
            // Apply crossover and mutation operations to create new individuals
            for (int i = 1; i < populationSize ; i +=2 ) {
                int parent1Index = random.nextInt(populationSize);
                int parent2Index = random.nextInt(populationSize);
                boolean[][] child1 = new boolean[C][V];
                boolean[][] child2 = new boolean[C][V];
                // Crossover
                if (random.nextInt(100) < crossoverProb) {
                    for (int row = 0; row < C; row++) {
                        for (int col =0 ; col< V ; col++)
                        {
                            if (row < (C/2) )
                            {
                                child1[row][col] = population[parent1Index] [row] [col];
                            }
                            else 
                            {
                                child1[row][col] = population[parent2Index] [row] [col];
                            }
                        }
                    }
                    if (isValid(child1)) {
                        newPopulation.add(child1);
                    }
                    if (isValid(child2)) {
                        newPopulation.add(child2);
                    }
                    for (int row = 0; row < C; row++) {
                        for (int col =0 ; col< V ; col++)
                        {
                            if (col < (V/2) )
                            {
                                child1[row][col] = population[parent1Index] [row] [col];
                            }
                            else 
                            {
                                child1[row][col] = population[parent2Index] [row] [col];
                            }
                        }
                        if (isValid(child1)) {
                            newPopulation.add(child1);
                        }
                        if (isValid(child2)) {
                            newPopulation.add(child2);
                        }
                    }
                }
                // Mutation
                if (random.nextInt(100) < mutationProb) {

                    for (int row = 0; row <Math.min(100,C) ;row++) {
                        for (int col =0 ; col<Math.min(100,V); col++)
                        {
                            child1 = new boolean[C][V];
                            child2 = new boolean[C][V];
                            
                            // Perform deep copy for child1
                            for (int roww = 0; roww < C; roww++) {
                                for (int coll = 0; coll < V; coll++) {
                                    child1[roww][coll] = population[i][roww][coll];
                                }
                            }
                            
                            // Perform deep copy for child2
                            for (int roww = 0; roww < C; roww++) {
                                for (int coll = 0; coll < V; coll++) {
                                    child2[roww][coll] = population[i][roww][coll];
                                }
                            }
                            
                            child1[row][col] = !population[i][row][col];
                            if (isValid(child1)) {
                                newPopulation.add(child1);
                            }

                            child2[row][col] = ! population[i-1][row][col];
                            if (isValid(child2)) {
                                newPopulation.add(child2);
                            }
                        }
                    }
                }
           }
            double[] scores = new double[newPopulation.size()];
            for (int i = 0; i < newPopulation.size(); i++) {
                scores[i] = ScoreEstimate(newPopulation.get(i));
            }
            Arrays.sort(scores);
            double leastAvalibleScore = scores[scores.length - populationSize];
            // Create a new list to hold the selected individuals
            // Update newPopulation with the selected individuals
            int itr =0 ;
            boolean[][][] nextGeneration = new boolean[populationSize][C][V];
            for (boolean[][] individual : newPopulation) {
                
                    double score = ScoreEstimate(individual);
                    if (score > leastAvalibleScore) {
                        nextGeneration[itr] = individual;
                        itr++;
                    }

                }
                for (boolean[][] individual : newPopulation) {
                
                    double score = ScoreEstimate(individual);
                    if (score == leastAvalibleScore) {
                        nextGeneration[itr] = individual;
                        itr++;
                    }
                    if (itr == populationSize)
                    {
                        break;
                    }
                }
            // Update the population with the next generation
            population = nextGeneration;
            double mx  = scores[scores.length-1];
            // Check if the maximum fitness score threshold is reached
            double IncreasingRatio = ( mx/(oldmax) )-1;
            if (IncreasingRatio < minFitnessScoreRatioThreshold) {
                break;
            }
            oldmax=mx;
        }
     
        return population;
    }

    public static boolean[][] hillClimbing() {
        boolean[][] Sol = new boolean [C] [V];
        for (int i=0 ; i<C ;i++){
            for (int j=0 ; j<V ;j++)
            {
                Sol[i][j]=false;
            }
        }
        double currentScore = ScoreEstimate(Sol);
        int maxIter = 100000; // maximum iterations
        for (int i = 0; i < maxIter; i++) {
            boolean[][] nextsol = Sol.clone();
            int row = new Random().nextInt(C);
            int col = new Random().nextInt(V);
            nextsol[row][col] = !nextsol[row][col];
            if (isValid(nextsol)) {
                double nextScore = ScoreEstimate(nextsol);
                if (nextScore > currentScore) {
                    Sol = nextsol;
                    currentScore = nextScore;
                }
            }
        }

        return Sol;
    }

    // This function here gives a score to the solution 2-D array which is provided as a parameter
    public static double Scoringfunction ( boolean [] [] sol  ) {
        double totalScore = 0;
        int totalRequests = 0;
    
        // Iterate through each request
        for (int endpointId = 0; endpointId < requestFrequency.length; endpointId++) {
            for (int vidId = 0; vidId < requestFrequency[endpointId].length; vidId++) {
                int currentRequestFrequency = requestFrequency[endpointId][vidId];
                totalRequests += currentRequestFrequency;
    
                int bestLatency = dataCenterCost[endpointId]; // Initialize with data center cost
    
                // Iterate through each cache
                for (int cacheId = 0; cacheId < sol.length; cacheId++) {
                    if (sol[cacheId][vidId]) {
                        int cacheLatency = serverCost[endpointId][cacheId];
                        if (cacheLatency < bestLatency) {
                            bestLatency = cacheLatency;
                        }
                    }
                }
    
                totalScore += (currentRequestFrequency * (dataCenterCost[endpointId] - bestLatency));
            }
        }
    
        double scorePerRequest = totalScore / totalRequests;
        double score = scorePerRequest * 1000;
    
        return score;
    }
    public static double ScoreEstimate (boolean [][] sol )
    {
        int score = 0 ;
        for (int i=0 ; i<C; i++)
        {
            for (int j = 0 ; j<V ; j++)
            {
                if (sol [i][j])
                {
                    score += CellsScore [i][j];
                }
            }
        }
        return score/R;
    }
}
