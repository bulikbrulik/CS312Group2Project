using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DPTSP
{
    class Memo : IEquatable<Memo>
    {
        public Subset subset;
        public int end;

        public Memo(Subset subset, int end)
        {
            this.subset = subset;
            this.end = end;
        }

        public bool Equals(Memo other)
        {
            if(ReferenceEquals(null, other))
            {
                return false;
            }
            return subset == other.subset && end == other.end;
        }

        public override bool Equals(object obj)
        {
            return this.Equals(obj as Memo);
        }

        public static bool operator ==(Memo left, Memo right)
        {
            if(ReferenceEquals(null, left))
            {
                return ReferenceEquals(null, right);
            }
            return left.Equals(right);
        }

        public static bool operator !=(Memo left, Memo right)
        {
            return !(left == right);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                if(!ReferenceEquals(null, subset))
                {
                    hash = 23 * subset.GetHashCode();
                }
                hash = 23 * end;
                return hash;
            }
        }
    }
}
