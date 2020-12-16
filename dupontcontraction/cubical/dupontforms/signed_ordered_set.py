"""
Class for ordered sets, returns the shuffle sign whenever adding elements
(added from the right) or joining sets together.
"""

class SignedOrderedSet():
    
    def __init__(self, iterable=None):
        self.start = None  # smallest element
        self.end = None  # biggest elements
        self.set = dict()  # key => [previous, next]
        self._sign = 1  # permutation sign
        
        if iterable is not None:
            for el in iterable:
                self.add_el(el)
    
    
    def __len__(self):
        return len(self.set)
    
    
    def add_el(self, el):
        if len(self) == 0:
            self.start, self.end = el, el
            self.set[el] = [None, None]
            self._sign *= 1
            return
        
        if el in self.set:
            self._sign = 0  # repeated elements give 0
            return
        
        if el > self.end:
            self.set[self.end][1] = el
            self.set[el] = [self.end, None]
            self.end = el
            self._sign *= 1
            return
        
        if el < self.start:
            sign = (-1)**len(self)
            self.set[self.start][0] = el
            self.set[el] = [None, self.start]
            self.start = el
            self._sign *= sign
            return
        
        sign = 1
        before, after = None, self.end
        while after > el:
            sign *= -1
            before = after
            after = self.set[after][0]
        self.set[el] = [after, before]
        self.set[after][1], self.set[before][0] = el, el
        self._sign *= sign
        return
    
    
    def __iter__(self):
        curr = self.start
        while curr != self.end:
            yield curr
            curr = self.set[curr][1]
        yield curr
        return
    
    
    def __add__(self, other):
        out = SignedOrderedSet(self.set)
        for el in other:
            out.add_el(el)
        return out
    
    
    def sign(self):
        """
        Output the sign and reset sign to 1.
        """
        out = self._sign
        self._sign = 1
        return out
