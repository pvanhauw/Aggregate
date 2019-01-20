'''
Created on Jan 19, 2019

@author: pierre
'''

# https://docs.python.org/2/library/operator.html 

import math
assert(math.sqrt)

class Point(object):
    default_init_float = 0. 
    def __init__(self, x=default_init_float , y=default_init_float, z=default_init_float ):
        self.x = x 
        self.y = y 
        self.z = z 
        
    def __str__(self):
        return '%f %f %f'% (self.x, self.y,self.z)
        
    def __add__(self, point):
        return Point(point.x +self.x, point.y +self.y, point.z +self.z  )
    
    def __iadd__(self, point):
        self = self + point 
        return self
    
    def __sub__(self, point):
        return Point(point.x - self.x, point.y - self.y, point.z - self.z  )
    
    def __isub__(self, point):
        self = self - point 
        return self
    
    def __mul__(self, scalar ):
        return Point(self.x * scalar, self.y * scalar , self.z * scalar )
    
    def __imul__(self, scalar):
        self = self * scalar
        return self

    def __truediv__(self, scalar ):
        return Point(self.x / scalar, self.y / scalar , self.z / scalar )
    
    def __itruediv__(self, scalar):
        self = self / scalar
        return self
    
    # cross prod  -> " ^ " 
    def __xor__(self, p):
        x = self.y * p.z - self.z * p.y  
        y = self.z * p.x - self.x * p.z  
        z = self.x * p.y - self.y * p.x  
        return Point(x,y,z ) 
    
    def Normalize(self):
        tmp = self.x * self.x + self.y * self.y + self.z * self.z
        try :
            norm = math.sqrt(tmp)
            # self /= norm  # doesnot work 
            self.x = self.x / norm 
            self.y = self.y / norm 
            self.z = self.z / norm 
        except:
            print("div per zero for point "+str(self))
        return self 
    
    def Norm(self):
        tmp = self.x * self.x + self.y * self.y + self.z * self.z
        return math.sqrt(tmp)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    