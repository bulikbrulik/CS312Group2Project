using System;
using System.Collections;
using System.Linq;
using System.Text;

namespace TSP
{
    class TSPState
    {
        private ArrayList path;
        private double lowerBound;
        private double priority;
        private double[,] costMatrix;

        public TSPState(ref ArrayList newPath, ref double newLowerBound, ref double[,] newCostMatrix)
        {
            path = newPath;
            lowerBound = newLowerBound;
            costMatrix = newCostMatrix;
            priority = double.MaxValue;
        }

        public ArrayList getPath() { return path; }
        public void addCityToPath(City city) { path.Add(city); }
        public double getPriority() { return priority; }
        public void setPriority(double prio) { priority = prio; }
        public double getLowerBound() { return lowerBound; }
        public void setLowerBound(double bound) { lowerBound = bound; }
        public double[,] getCostMatrix() { return costMatrix; }
    }
}
