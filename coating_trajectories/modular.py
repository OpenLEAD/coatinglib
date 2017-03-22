class Z:
    mod = 5
    def __init__(self,n):
            self.n = int(n) % Z.mod
    def __add__(self,other):
            return Z(self.n + other)
    def __mul__(self,other):
            return Z(self.n*other)
    def __div__(self,other):
            return self * Z(other).inv()
    def __sub__(self,other):
            return Z(self.n - other)
    def __int__(self):
            return self.n % Z.mod
    def __radd__(self,other):
            return Z(self.n + other)
    def __rsub__(self,other):
            return Z(other - self.n)
    def __rmul__(self,other):
            return Z(self.n*other)
    def __rdiv__(self,other):
            return self.inv() * other
    def __repr__(self):
            return str(self.n % Z.mod)
    def inv(self):
        t, newt = 0, 1
        r, newr = Z.mod, self.n

        while newr != 0:
            quotient= r / newr
            t, newt = newt, t - quotient * newt
            r, newr = newr, r - quotient * newr

        if r > 1:
            raise ValueError(str(self) + " has no inverse on Z/"+str(Z.mod)+"Z.")

        return Z(t)
