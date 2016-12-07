using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;
using DPTSP;
using System.Linq;

namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }
        //************************************************Helper Functions**********************************************************************************
        //Creates the initial cost matrix by looping through the matrix and the diagonal to inifinity and the costs to the edge costs. Takes O(n^2) time and space
        private TSPState createInitialState()
        {
            //Loop throught the matrix setting all the intial values. Takes O(n^2) space and time
            double[,] matrix = new double[Cities.Length, Cities.Length];
            for (int i = 0; i < Cities.Length; i++)
            {
                for(int j = 0; j < Cities.Length; j++)
                {
                    if(i == j)
                    {
                        matrix[i,j] = double.MaxValue;
                    }
                    else
                    {
                        matrix[i, j] = Cities[i].costToGetTo(Cities[j]);
                    }
                }
            }
            //Set the path to just the starting node. Takes O(1) time and space
            ArrayList path = new ArrayList();
            path.Add(Cities[0]);

            //Calculate the lower bound. Takes O(n^2) time and O(1) space
            double lowerBound = reduceMatrix(ref matrix);

            return new TSPState(ref path, ref lowerBound, ref matrix);

        }

        //Goes through each row and column in the matrix and subtracts the lowest value of that row from the other entries in the row. 
        //Takes O(n^2) since it has to go through every index in the matrix two times so it is technically O(4n^2).
        private double reduceMatrix(ref double[,] matrix)
        {
            double lowerBound = 0;
            //Go through each row to reduce
            for (int r = 0; r < Cities.Length; r++)
            {
                double minVal = double.MaxValue;
                for(int c = 0; c < Cities.Length; c++)
                {
                    if(matrix[r,c] < minVal)
                    {
                        minVal = matrix[r, c];
                    }
                }

                if(minVal != 0 && minVal != double.MaxValue)
                {
                    lowerBound += minVal;
                    for(int c = 0; c < Cities.Length; c++)
                    {
                        if(matrix[r,c] != double.MaxValue)
                        {
                            matrix[r, c] -= minVal;
                        }
                    }
                }
            }
            //Go through each column to reduce
            for (int c = 0; c < Cities.Length; c++)
            {
                double minVal = double.MaxValue;
                for (int r = 0; r < Cities.Length; r++)
                {
                    if (matrix[r, c] < minVal)
                    {
                        minVal = matrix[r, c];
                    }
                }

                if (minVal != 0 && minVal != double.MaxValue)
                {
                    lowerBound += minVal;
                    for (int r = 0; r < Cities.Length; r++)
                    {
                        if (matrix[r, c] != double.MaxValue)
                        {
                            matrix[r, c] -= minVal;
                        }
                    }
                }
            }


            return lowerBound;
        }
        
        //Creates a key for the states to use in the priority queue. O(1) for space and time since it only performs addition and multiplication and no extra space is created
        private double calculateKey(int citiesLeft, double lowerBound)
        {
            if(citiesLeft < 1) { return lowerBound; }
            else { return lowerBound / (Cities.Length - citiesLeft); }
        }

        //Creates an initial path greedily. Time complexity is O(n^2) because it goes through each city (or edge) in the list to find the smallest edge cost.
        //Space complexity is O(n) because it creates a list the size of the number of cities.
        private double createGreedyBssf()
        {
            Route = new ArrayList();
            Route.Add(Cities[0]);
            int currCityIndex = 0;

            while(Route.Count < Cities.Length)
            {
                double minVal = double.MaxValue;
                int minIndex = 0;

                for(int i = 0; i < Cities.Length; i++)
                {
                    if(currCityIndex != i)
                    {
                        if (!Route.Contains(Cities[i]))
                        {
                            double currCost = Cities[currCityIndex].costToGetTo(Cities[i]);
                            if(currCost < minVal)
                            {
                                if(Route.Count == Cities.Length-1 && Cities[i].costToGetTo(Cities[0]) == double.MaxValue)
                                {
                                    continue;
                                }
                                minVal = currCost;
                                minIndex = i;
                            }
                        }
                    }
                }
                currCityIndex = minIndex;
                Route.Add(Cities[currCityIndex]);
            }
            bssf = new TSPSolution(Route);
            return bssf.costOfRoute();
        }
        
        /**
       * Helper Function to initially set up a cost matrix at a current state
       * Time Complexity: O(n) as it iterates over one row and one column in the matrix and n would be the length of 
       * the row/column which is the number of cities in the graph, or rather |V|
       * Space Complexity: O(1) as all the data is passed by reference and the function does not create extra
       * data structures that depend on the size of the input.
       */
        void setUpMatrix(ref double[,] costMatrix, int indexOfParent, int indexOfChild, ref double lowerBound)
        {
            if (costMatrix[indexOfParent, indexOfChild] != double.MaxValue)
                lowerBound += costMatrix[indexOfParent, indexOfChild];
            // Make sure to set all costs coming from the currState to infinity
            for (int column = 0; column < Cities.Length; column++)
            {
                costMatrix[indexOfParent, column] = double.MaxValue;
            }
            // Make sure to set all costs coming into the child State to infinity
            for (int row = 0; row < Cities.Length; row++)
            {
                costMatrix[row, indexOfChild] = double.MaxValue;
            }
            // Make sure to set the cost of going from child state back to parent to infinity as we don't want cycles
            costMatrix[indexOfChild, indexOfParent] = double.MaxValue;
        }

        //**********************************************BBAlgorithm*********************************************************************************
        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for a branch and bound solver here.
            //Initalize variables. Takes O(1) space and time
            int numOfCitiesLeft = Cities.Length;
            int numOfSolutions = 0;
            int numOfStatesCreated = 0;
            int numOfStatesNotExpanded = 0;
            
            //Initalize time variable for stopping the algorithm after the default of 60 seconds. Takes O(1) space and time
            DateTime start = DateTime.Now;
            DateTime end = start.AddSeconds(time_limit / 1000);

            //Create initial root state and set its priority to its lower bound. Takes O(n^2) space and time as discussed above
            TSPState initialState = createInitialState();
            numOfStatesCreated++;
            initialState.setPriority(calculateKey(numOfCitiesLeft - 1, initialState.getLowerBound()));

            //Create initial BSSF greedily
            double bssfBound = createGreedyBssf();

            PriorityQueue queue = new PriorityQueue(Cities.Length);
            queue.insert(initialState);

            // Branch and Bound until the queue is empty, we have exceeded the time limit, or we found the optimal solution
            /* This loop will have a iterate 2^n times approximately with expanding and pruning for each state, then for each state it
            does O(n^2) work by reducing the matrix, so over all O((n^2)*(2^n)) time and space as well as it creates a nxn 
            matrix for each state*/
            while (!queue.isEmpty() && DateTime.Now < end && queue.getMinLB() != bssfBound)
            {
                // Grab the next state in the queue
                TSPState currState = queue.deleteMin();

                // check if lower bound is less than the BSSF, else prune it
                if (currState.getLowerBound() < bssfBound)
                {
                    // Branch and create the child states
                    for (int i = 0; i < Cities.Length; i++)
                    {
                        // First check that we haven't exceeded the time limit
                        if (DateTime.Now >= end)
                            break;

                        // Make sure we are only checking cities that we haven't checked already
                        if (currState.getPath().Contains(Cities[i]))
                            continue;

                        // Create the State
                        double[,] oldCostMatrix = currState.getCostMatrix();
                        double[,] newCostMatrix = new double[Cities.Length, Cities.Length];
                        // Copy the old array in the new one to modify the new without affecting the old
                        for (int k = 0; k < Cities.Length; k++)
                        {
                            for (int l = 0; l < Cities.Length; l++)
                            {
                                newCostMatrix[k, l] = oldCostMatrix[k, l];
                            }
                        }
                        City lastCityinCurrState = (City)currState.getPath()[currState.getPath().Count - 1];
                        double oldLB = currState.getLowerBound();
                        setUpMatrix(ref newCostMatrix, Array.IndexOf(Cities, lastCityinCurrState), i, ref oldLB);
                        double newLB = oldLB + reduceMatrix(ref newCostMatrix);
                        ArrayList oldPath = currState.getPath();
                        ArrayList newPath = new ArrayList();
                        foreach (City c in oldPath)
                        {
                            newPath.Add(c);
                        }
                        newPath.Add(Cities[i]);
                        TSPState childState = new TSPState(ref newPath, ref newLB, ref newCostMatrix);
                        numOfStatesCreated++;

                        // Prune States larger than the BSSF
                        if (childState.getLowerBound() < bssfBound)
                        {
                            City firstCity = (City)childState.getPath()[0];
                            City lastCity = (City)childState.getPath()[childState.getPath().Count - 1];
                            double costToLoopBack = lastCity.costToGetTo(firstCity);

                            // If we found a solution and it goes back from last city to first city
                            if (childState.getPath().Count == Cities.Length && costToLoopBack != double.MaxValue)
                            {
                                childState.setLowerBound(childState.getLowerBound() + costToLoopBack);
                                bssf = new TSPSolution(childState.getPath());
                                bssfBound = bssf.costOfRoute();
                                numOfSolutions++;
                                numOfStatesNotExpanded++; // this state is not expanded because it is not put on the queue
                            }
                            else
                            {
                                // Set the priority for the state and add the new state to the queue
                                numOfCitiesLeft = Cities.Length - childState.getPath().Count;
                                childState.setPriority(calculateKey(numOfCitiesLeft, childState.getLowerBound()));
                                queue.insert(childState);
                            }
                        }
                        else
                        {
                            numOfStatesNotExpanded++; // States that are pruned are not expanded
                        }
                    }
                }
                currState = null;
            }
            numOfStatesNotExpanded += queue.getSize(); // if the code terminated before queue is empty, then those states never got expanded
            Console.WriteLine("Number of states generated: " + numOfStatesCreated);
            Console.WriteLine("Number of states not Expanded: " + numOfStatesNotExpanded);
            Console.WriteLine("Max Number of states put in queue: " + queue.getMaxNumOfItems());
            end = DateTime.Now;
            TimeSpan diff = end - start;
            double seconds = diff.TotalSeconds;
            results[COST] = System.Convert.ToString(bssf.costOfRoute());    // load results into array here, replacing these dummy values
            results[TIME] = System.Convert.ToString(seconds);
            results[COUNT] = System.Convert.ToString(numOfSolutions);

            return results;
        }

       
        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////
        //************************************************************************Greedy Algorithm*****************************************************************************************************
        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            //Initialize the timer
            DateTime start = DateTime.Now;
            DateTime end = start.AddSeconds(time_limit / 1000);
            Random random = new Random();
            Route = new ArrayList();

            results[COST] = "not implemented";    // load results into array here, replacing the default values
            while (DateTime.Now < end && Route.Count < Cities.Length)
            {
                // Create variables to track progress
                Route.Add(Cities[random.Next(0, Cities.Length)]);
                int currCityIndex = 0;

                // While we haven't added |V| edges to our route
                while (Route.Count < Cities.Length)
                {
                    double minValue = double.MaxValue;
                    int minIndex = 0;

                    // Loop over all the cities and find the one with min cost to get to
                    for (int i = 0; i < Cities.Length; i++)
                    {
                        // We don't want to be checking ourselves or ones that are already in the route because it won't be a tour
                        if (currCityIndex != i && !Route.Contains(Cities[i]))
                        {
                            double tempValue = Cities[currCityIndex].costToGetTo(Cities[i]);
                            if (tempValue < minValue)
                            {
                                if (Route.Count == Cities.Length - 1 && Cities[i].costToGetTo(Cities[0]) == double.MaxValue)
                                {
                                    continue;
                                }
                                minValue = tempValue;
                                minIndex = i;
                            }
                        }
                    }

                    // Add the min edge to the Route by adding the destination city
                    currCityIndex = minIndex;
                    Route.Add(Cities[currCityIndex]);
                }

                City lastCity = (City)(Route[Route.Count - 1]);
                City firstCity = (City)(Route[0]);
                if (lastCity.costToGetTo(firstCity) == double.MaxValue)
                {
                    Route.Clear();
                }

            }

            int numOfSolutions = 0;
            // Display the results
            bssf = new TSPSolution(Route);
            end = DateTime.Now;
            TimeSpan diff = end - start;
            double seconds = diff.TotalSeconds;
            results[COST] = System.Convert.ToString(bssf.costOfRoute());    // load results into array here, replacing these dummy values
            results[TIME] = System.Convert.ToString(seconds);
            results[COUNT] = System.Convert.ToString(numOfSolutions);

            return results;
        }

        // good speed, bad memory
        public string[] fancySolveProblem()
        {
            string[] results = new string[3];
            Stopwatch timer = new Stopwatch();
            int count = 0;
            timer.Start();

            double[][] costs = getCostMatrix();
            int n = Cities.Length;
            int npow = (int)Math.Pow(2, n);
            double[][] C = new double[npow][];
            C[1] = new double[n];
            C[1][0] = 0;
            for(int s = 2; s <= n; s++) // go through all subset sizes
            {
                for(int S = 3; S < npow; S+=2) // go through all subsets with city 0
                {
                    if(numberOfSetBits(S) == s) // |S| = s
                    {
                        if(C[S] == null)
                        {
                            C[S] = new double[n];
                        }
                        C[S][0] = double.MaxValue; // can't go this way
                        for(int j = 1; j < n; j++) // go through each city (exclude 0)
                        {
                            int jbit = (int)Math.Pow(2, j); // get the bit for the city
                            if((S & jbit) == jbit) // if this city is in the subset
                            {
                                double min = double.MaxValue;
                                for(int i = 0; i < n; i++) // go through each city (exclude 0)
                                {
                                    if(i != j)
                                    {
                                        int ibit = (int)Math.Pow(2, i); // get the bit for the city
                                        if((S & ibit) == ibit) // if this city is in the subset
                                        {
                                            int withoutj = (S & (~jbit));
                                            double cost = C[withoutj][i] + costs[i][j];
                                            if(cost < 0)
                                            {
                                                cost = double.MaxValue;
                                            }
                                            min = (cost < min) ? cost : min;
                                        }
                                    }
                                }
                                C[S][j] = min;
                            }
                        }
                    }
                }
            }

            int tour = npow - 1;
            double totalCost = double.MaxValue;
            ArrayList path = new ArrayList();
            int end = 0;
            while (tour != 1)
            {
                double min = double.MaxValue;
                int newEnd = n + 1;
                for(int j = 1; j < n; j++)
                {
                    int jbit = (int)Math.Pow(2, j);
                    if((tour & jbit) == jbit)
                    {

                        double cost = C[tour][j] + costs[j][end];
                        if(cost < min)
                        {
                            min = cost;
                            newEnd = j;
                            if(tour == npow - 1)
                            {
                                totalCost = cost;
                            }
                        }

                    }
                }
                end = newEnd;
                int endBit = (int)Math.Pow(2, end);
                tour = (tour & (~endBit));
                path.Add(Cities[end]);
            }
            path.Add(Cities[0]);
            path.Reverse();
            count = 1;

            timer.Stop();
            bssf = new TSPSolution(path);
            results[COST] = totalCost.ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();
            return results;
        }

        int numberOfSetBits(int i)
        {
            i = i - ((i >> 1) & 0x55555555);
            i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
            return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
        }

        // bad speed, good memory
        /*public string[] fancySolveProblem()
        {
            int count = 0;
            string[] results = new string[3];
            defaultSolveProblem();
            double[][] costs = getCostMatrix();

            Stopwatch timer = new Stopwatch();
            timer.Start();

            Dictionary<Memo, double> memos = new Dictionary<Memo, double>();
            List<Subset> subsets = GetSubsets();

            foreach(Subset subset in subsets)
            {
                Memo baseCase = new Memo(subset, 0);
                if(subset.count() == 1)
                {
                    memos[baseCase] = 0;
                }
                else
                {
                    memos[baseCase] = double.MaxValue;
                }
            }

            for(int m = 2; m <= Cities.Length; m++)
            {
                foreach(Subset subset in subsets)
                {
                    if(subset.count() == m)
                    {
                        List<int> S = subset.getCities();
                        foreach(int j in S)
                        {
                            if(j != 0)
                            {
                                Memo current = new Memo(subset, j);
                                double min = double.MaxValue;
                                foreach(int k in S)
                                {
                                    if(k != j)
                                    {
                                        Memo sub = new Memo(subset.minus(j), k);
                                        double cost = memos[sub] + costs[k][j];
                                        min = (cost < min) ? cost : min;
                                    }
                                }
                                memos[current] = min;                       
                            }
                        }
                    }
                }
            }

            Subset tour = new Subset(new List<int>());
            for(int i = 0; i < Cities.Length; i++)
            {
                tour.add(i);
            }
            double totalCost = double.MaxValue;
            Memo from = new Memo(new Subset(new List<int>(new int[] {int.MaxValue})), int.MaxValue);
            ArrayList path = new ArrayList();
            int end = 0;
            bool stop = false;
            while (tour != new Subset(new List<int>(new int[] {0})) && !stop)
            {
                double min = double.MaxValue;
                List<int> cities = tour.getCities();
                foreach (int j in cities)
                {
                    if(j != 0)
                    {
                        Memo temp = new Memo(tour, j);
                        if(tour.count() == Cities.Length)
                        {
                            double cost = memos[temp] + costs[j][end];
                            if(cost < min)
                            {
                                totalCost = cost;
                                min = cost;
                                from = temp;
                                //stop = true;
                            }
                        }
                        else
                        {
                            double cost = memos[temp] + costs[j][end];
                            if(cost < min)
                            {
                                min = cost;
                                from = temp;
                            }
                        }
                    }
                }
                tour = tour.minus(from.end);
                path.Add(Cities[from.end]);
                end = from.end;
                Console.WriteLine(from.end);
            }
            path.Add(Cities[0]);
            path.Reverse();
            count = 1;

            timer.Stop();

            bssf = new TSPSolution(path);
            results[COST] = totalCost.ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        public List<Subset> GetSubsets()
        {
            List<Subset> subsets = new List<Subset>();
            Subset start = new Subset(new List<int>(new int[] {0}));
            subsets.Add(start);
            bool changed = true;
            while(changed)
            {
                changed = false;
                for(int i = 0; i < subsets.Count; i++)
                {
                    Subset subset = subsets[i];
                    for(int j = 0; j < Cities.Length; j++)
                    {
                        if(!subset.contains(j))
                        {
                            Subset bigger = new Subset(subset);
                            bigger.add(j);
                            if(!subsets.Contains(bigger))
                            {
                                subsets.Add(bigger);
                                changed = true;
                            }
                        }
                    }
                }
            }
            return subsets;
        }*/

        public double[][] getCostMatrix()
        {
            double[][] costs = new double[Cities.Length][];
            for(int from = 0; from < Cities.Length; from++)
            {
                costs[from] = new double[Cities.Length];
                for(int to = 0; to < Cities.Length; to++)
                {
                    costs[from][to] = (from == to) ? double.MaxValue : Cities[from].costToGetTo(Cities[to]);
                }
            }
            return costs;
        }
        #endregion
    }

}
