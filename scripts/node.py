class Node:
    def __init__(self, name, pop_size=None, left=None, right=None):
        self.name = name
        self.pop_size = pop_size
        self.left = left 
        self.right = right
        self.parent = None
        self.marginal = None
        self.transition = None # forward transition probability matrix; P(node | parent)

        if self.left:
            self.left.parent = self

        if self.right:
            self.right.parent = self

    def children(self):
        if not self.is_parent():
            return([self.name])
        else:
            child_list = []

            if self.left:
                child_list.extend(self.left.children())

            if self.right:
                child_list.extend(self.right.children())

            return(child_list)
    
    def is_parent(self):
        return(self.left is not None or self.right is not None)

    def is_ancestor(self):
        return(self.parent is None)