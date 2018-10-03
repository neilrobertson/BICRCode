# Originally from: https://bitbucket.org/james_taylor/bx-python/src/d9c88c9359a0/lib/bx/intervals/intersection.pyx

import math
import sys
import operator
import random

"""
Data structure for performing intersect queries on a set of intervals which
preserves all information about the intervals (unlike bitset projection methods).

:Authors: James Taylor (james@jamestaylor.org),
          Ian Schenk (ian.schenck@gmail.com),
          Brent Pederson (bpederse@gmail.com),
          Tony McBryan (tony@mcbryan.co.uk)
"""

# Historical note:
#    This module original contained an implementation based on sorted endpoints
#    and a binary search, using an idea from Scott Schwartz and Piotr Berman.
#    Later an interval tree implementation was implemented by Ian for Galaxy's
#    join tool (see `bx.intervals.operations.quicksect.py`). This was then
#    converted to Cython by Brent, who also added support for
#    upstream/downstream/neighbor queries. This was modified by James to
#    handle half-open intervals strictly, to maintain sort order, and to
#    implement the same interface as the origianl Intersecter. This was 
#    modified by Tony to convert back to pure python due to quadratic behaviour
#    under python 2.7.4 (that did not exist under 2.7.1) that only occurred under
#    the Cython version.

class EmptyNode:
    def __init__(self):
        self.croot = self

EmptyNode = EmptyNode()

RAND_MAX = sys.maxint
nlog = -1.0 / math.log(0.5)


class IntervalNode:
    """
    A single node of an `IntervalTree`.
    
    NOTE: Unless you really know what you are doing, you probably should us
          `IntervalTree` rather than using this directly. 
    """


    
    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __init__(self, start, end, interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority   = math.ceil(nlog * math.log(-1.0/(1.0 * random.uniform(0,RAND_MAX-1)/RAND_MAX - 1)))
        self.start      = start
        self.end       = end
        self.interval   = interval
        self.maxend    = end
        self.minstart   = start
        self.minend    = end
        self.cleft       = EmptyNode
        self.cright      = EmptyNode
        self.croot       = EmptyNode
        
    def insert(self, start, end, interval):
        """
        Insert a new IntervalNode into the tree of which this node is
        currently the root. The return value is the new root of the tree (which
        may or may not be this node!)
        """
        croot = self
        # If starts are the same, decide which to add interval to based on
        # end, thus maintaining sortedness relative to start/end
        decision_endpoint = start
        if start == self.start:
            decision_endpoint = end
        
        if decision_endpoint > self.start:
            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright.insert( start, end, interval )
            else:
                self.cright = IntervalNode( start, end, interval )
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()
        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft.insert( start, end, interval)
            else:
                self.cleft = IntervalNode( start, end, interval)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()
    
        croot.set_ends()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    def rotate_right(self):
        croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_ends()
        return croot

    def rotate_left(self):
        croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_ends()
        return croot

    def set_ends(self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode: 
            self.maxend = max(self.end, self.cright.maxend, self.cleft.maxend)
            self.minend = min(self.end, self.cright.minend, self.cleft.minend)
            self.minstart = min(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxend = max(self.end, self.cright.maxend)
            self.minend = min(self.end, self.cright.minend)
            self.minstart = min(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxend = max(self.end, self.cleft.maxend)
            self.minend = min(self.end, self.cleft.minend)
            self.minstart = min(self.start, self.cleft.minstart)
        
    def intersect(self, start, end):
        """
        given a start and a end, return a list of features
        falling within that range
        """
        results = []
        self._intersect( start, end, results )
        return results

    find = intersect
    
    def _intersect(self, start, end, results ):
        # Left subtree
        if self.cleft is not EmptyNode and self.cleft.maxend > start:
            self.cleft._intersect( start, end, results )
        # This interval
        if ( self.end > start ) and ( self.start < end ):
            results.append( self.interval )
        # Right subtree
        if self.cright is not EmptyNode and self.start < end:
            self.cright._intersect( start, end, results )
    
    def _seek_left(self, position, results, n, max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend + max_dist < position:
            return
        if self.minstart > position:
            return

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cright is not EmptyNode:
            self.cright._seek_left(position, results, n, max_dist)

        if -1 < position - self.end < max_dist:
            results.append(self.interval)

        # TODO: can these conditionals be more stringent?
        if self.cleft is not EmptyNode:
                self.cleft._seek_left(position, results, n, max_dist)


    
    def _seek_right(self, position, results, n, max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxend, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode: 
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    
    def left(self, position, n=1, max_dist=2500):
        """
        find n features with a start > than `position`
        f: a Interval object (or anything with an `end` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left( position - 1, results, n, max_dist )
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        return r[:n]

    def right(self, position, n=1, max_dist=2500):
        """
        find n features with a end < than position
        f: a Interval object (or anything with a `start` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        results = []
        # use end + 1 becuase .right() assumes strictly right-of
        self._seek_right(position + 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        return r[:n]

    def traverse(self, func):
        self._traverse(func)

    def _traverse(self, func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self)
        if self.cright is not EmptyNode: self.cright._traverse(func)

    def get_interval(self):
        return self.interval

## ---- Wrappers that retain the old interface -------------------------------

class Interval:
    """
    Basic feature, with required integer start and end properties.
    a name, and any arbitrary data is sent in on the info keyword argument

    >>> from bx.intervals.intersection import Interval

    >>> f1 = Interval(23, 36)
    >>> f2 = Interval(34, 48, value={'chr':12, 'anno':'transposon'})
    >>> f2
    Interval(34, 48, value={'anno': 'transposon', 'chr': 12})

    """

    def __init__(self, start, end, value=None):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.end   = end      
        self.value = value

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        fstr += ")"
        return fstr

    def __cmp__(self, other):
        return cmp( self.start, other.start ) or cmp( self.end, other.end )

class IntervalTree:
    """
    Data structure for performing window intersect queries on a set of 
    of possibly overlapping 1d intervals.
    
    Usage
    =====
    
    Create an empry IntervalTree
    
    >>> from bx.intervals.intersection import Interval, IntervalTree
    >>> intersecter = IntervalTree()
    
    An interval is a start and end position and a value (possibly None).
    You can add any object as an interval:
    
    >>> intersecter.insert( 0, 10, "food" )
    >>> intersecter.insert( 3, 7, dict(foo='bar') )
    
    >>> intersecter.find( 2, 5 )
    ['food', {'foo': 'bar'}]
    
    If the object has start and end attributes (like the Interval class) there
    is are some shortcuts:
    
    >>> intersecter = IntervalTree()
    >>> intersecter.insert_interval( Interval( 0, 10 ) )
    >>> intersecter.insert_interval( Interval( 3, 7 ) )
    >>> intersecter.insert_interval( Interval( 3, 40 ) )
    >>> intersecter.insert_interval( Interval( 13, 50 ) )
    
    >>> intersecter.find( 30, 50 )
    [Interval(3, 40), Interval(13, 50)]
    >>> intersecter.find( 100, 200 )
    []
    
    Before/after for intervals
    
    >>> intersecter.before_interval( Interval( 10, 20 ) )
    [Interval(3, 7)]
    >>> intersecter.before_interval( Interval( 5, 20 ) )
    []

    """
    
    def __init__( self ):
        self.root = None
    
    # ---- Position based interfaces -----------------------------------------
    
    def insert( self, start, end, value=None ):
        """
        Insert the interval [start,end) associated with value `value`.
        """
        if self.root:
            self.root = self.root.insert( start, end, value )
        else:
            self.root = IntervalNode( start, end, value )
        
    def find( self, start, end ):
        """
        Return a sorted list of all intervals overlapping (start,end) (i.e. exclusive).
        """
        if self.root:
            return self.root.find( start, end )
        else:
            return []
    
    def find_inclusive(self, start, end):
        """
        Return a sorted list of all intervals overlapping [start, end] (i.e. inclusive)
        """
        return self.find(start-1, end+1)
        
    def find_ucsc(self, start, end):
        """
        Return a sorted list of all intervals overlapping [start, end) (i.e. inclusive start, exclusive end)
        This is the UCSC format for regions.
        """
        # note that this is not as simple as doign find_inclusive(start,end-1) as this misses cases where the
        # they overlap where the start 
        
        intervals = self.find_inclusive(start, end)
        
        # turns out this is the simplest hacky solution to this although a very slight performance penalty
        
        for interval in intervals[:]: # shallow copy the list for the traversal
            # filter any which are only overlapping with end of either interval
            if interval.start == end or start == interval.end:
                # the start and end are in the same location which means they are back to back but not technically touching
                intervals.remove(interval)
        return intervals
    
    def before( self, position, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie before `position` and are no
        further than `max_dist` positions aways
        """
        if self.root:
            return self.root.left( position, num_intervals, max_dist )
        else:
            return []

    def after( self, position, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie after `position` and are no
        further than `max_dist` positions aways
        """
        if self.root:
            return self.root.right( position, num_intervals, max_dist )
        else:
            return []

    # ---- Interval-like object based interfaces -----------------------------

    def insert_interval( self, interval ):
        """
        Insert an "interval" like object (one with at least start and end
        attributes)
        """
        self.insert( interval.start, interval.end, interval )

    def before_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie completely before `interval`
        and are no further than `max_dist` positions aways
        """
        return self.root.left( interval.start, num_intervals, max_dist )

    def after_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie comletey after `interval` and
        are no further than `max_dist` positions aways
        """
        return self.root.right( interval.end, num_intervals, max_dist )
    
    # ---- Old 'Intersecter' interface ----------------------------------------

    def add( self, start, end, value=None ):
        """
        Synonym for `insert`.
        """
        self.insert( start, end, value )
    
    def add_interval( self, interval ):
        """
        Synonym for `insert_interval`.
        """
        self.insert( interval )
    
    # ---- Miscellaneous utility function added by pzs ------------------------

    def traverse( self, func ):
        """
        traverse the tree, applying func to each node
        """
        self.root.traverse(func)

    def isEmpty( self ):
        """
        returns True if the tree currently has no intervals inserted
        """
        return self.root == None

# For backward compatibility
Intersecter = IntervalTree
