using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TSP
{
    class PriorityQueue
    {
        private int capacity;
        private int count;
        private int maxNum;
        private TSPState[] states;
        public PriorityQueue(int numOfNodes)
        {
            states = new TSPState[1000000];
            capacity = numOfNodes;
            count = 0;
            maxNum = 0;
        }
        public bool isEmpty()
        {
            return count == 0;
        }

        public int getSize()
        {
            return count;
        }
        
        public double getMinLB()
        {
            return states[1].getLowerBound();
        }

        /**
        * This method returns the TSPState of the element with the minimum value and removes it from the queue. 
        * Time Complexity: O(log(|V|)) because removing a node is constant time as we have its position in
        * the queue, then to read just the heap we just bubble up the min value which takes as long as 
        * the depth of the tree which is log(|V|), where |V| is the number of nodes
        * Space Complexity: O(1) because we don't create any extra variables that vary with the size of the input.
        */
        public TSPState deleteMin()
        {
            // grab the node with min value which will be at the root
            TSPState minValue = states[1];
            states[1] = states[count];
            count--;
            // fix the heap
            int indexIterator = 1;
            while (indexIterator <= count)
            {
                // grab left child
                int smallerElementIndex = 2 * indexIterator;

                // if child does not exist, break
                if (smallerElementIndex > count)
                    break;

                // if right child exists and is of smaller value, pick it
                if (smallerElementIndex + 1 <= count && states[smallerElementIndex + 1].getPriority() < states[smallerElementIndex].getPriority())
                {
                    smallerElementIndex++;
                }

                if (states[indexIterator].getPriority() > states[smallerElementIndex].getPriority())
                {
                    // set the node's value to that of its smaller child
                    TSPState temp = states[smallerElementIndex];
                    states[smallerElementIndex] = states[indexIterator];
                    states[indexIterator] = temp;
                }

                indexIterator = smallerElementIndex;
            }
            // return the min value
            return minValue;
        }

        /**
        * This function returns the maximum number of items ever put in the queue
        */
        public int getMaxNumOfItems()
        {
            return maxNum;
        }
        /**
        * This method updates the nodes in the queue after inserting a new node
        * Time Complexity: O(log(|V|)) as reording the heap works by bubbling up the min value to the top
        * which takes as long as the depth of the tree which is log|V|.
        * Space Complexity: O(1) as it does not create any extra variables that vary with the size of the input.
        */
        public void insert(TSPState newState)
        {
            // update the count
            count++;
            states[count] = newState;
            if (count > maxNum)
                maxNum = count;

            // as long as its parent has a larger value and have not hit the root
            int indexIterator = count;
            while (indexIterator > 1 && states[indexIterator / 2].getPriority() > states[indexIterator].getPriority())
            {
                // swap the two nodes
                TSPState temp = states[indexIterator / 2];
                states[indexIterator / 2] = states[indexIterator];
                states[indexIterator] = temp;

                indexIterator /= 2;
            }
        }
    }
}
