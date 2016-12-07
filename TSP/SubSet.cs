using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DPTSP
{
    public class Subset : IEquatable<Subset>
    {
        private List<int> cities;

        public Subset(List<int> cities)
        {
            this.cities = cities;
        }

        public Subset(Subset other)
        {
            cities = new List<int>();
            for(int i = 0; i < other.count(); i++)
            {
                add(other.get(i));
            }
        }

        public void add(int city)
        {
            if(!cities.Contains(city))
            {
                cities.Add(city);
            }
        }

        public int get(int index)
        {
            return cities[index];
        }

        public bool contains(int city)
        {
            return cities.Contains(city);
        }

        public int count()
        {
            return cities.Count;
        }

        public void remove(int city)
        {
            cities.Remove(city);
        }
        public List<int> getCities()
        {
            return cities;
        }

        public Subset minus(int city)
        {
            Subset other = new Subset(this);
            other.remove(city);
            return other;
        }

        public override String ToString()
        {
            StringBuilder builder = new StringBuilder();
            for(int i = 0; i < count(); i++)
            {
                builder.Append(get(i));
                if(i != count() - 1)
                {
                    builder.Append(", ");
                }
            }
            return builder.ToString();
        }

        public bool Equals(Subset other)
        {
            if(ReferenceEquals(null, other))
            {
                return false;
            }
            for(int i = 0; i < other.count(); i++)
            {
                if(!contains(other.get(i)))
                {
                    return false;
                }
            }
            for(int i = 0; i < count(); i++)
            {
                if(!other.contains(get(i)))
                {
                    return false;
                }
            }
            return true;
        }

        public override bool Equals(object obj)
        {
            return this.Equals(obj as Subset);
        }

        public static bool operator ==(Subset left, Subset right)
        {
            if(ReferenceEquals(null, left))
            {
                return ReferenceEquals(null, right);
            }
            return left.Equals(right);
        }

        public static bool operator !=(Subset left, Subset right)
        {
            return !(left == right);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                if(!ReferenceEquals(null, cities))
                {
                    foreach(int city in cities)
                    {
                        hash = hash * 23 * city;
                    }
                }
                return hash;
            }
        }
    }
}
