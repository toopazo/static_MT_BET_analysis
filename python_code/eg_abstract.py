# Python program showing
# abstract base class work

from abc import ABC, abstractmethod


class Polygon(ABC):
    @abstractmethod
    def noofsides(self):
        pass

    def classname(self):
        print('Base class method: type(self) %s ' % type(self))


class Triangle(Polygon):
    def noofsides(self):
        print("I have 3 sides")


class Pentagon(Polygon):
    def noofsides(self):
        print("I have 5 sides")


class Hexagon(Polygon):
    def noofsides(self):
        print("I have 6 sides")


class Quadrilateral(Polygon):
    def noofsides(self):
        print("I have 4 sides")

    def classname(self):
        super(Quadrilateral, self).classname()
        # print(super().classname())
        print('Quadrilateral class method: type(self) %s ' % type(self))


if __name__ == "__main__":
    # Driver code
    obj = Triangle()
    obj.noofsides()
    obj.classname()
    print('--')

    obj = Quadrilateral()
    obj.noofsides()
    obj.classname()
    print('--')

    obj = Pentagon()
    obj.noofsides()
    obj.classname()
    print('--')

    obj = Hexagon()
    obj.noofsides()
    obj.classname()
    print('--')

    assert(isinstance(obj, Polygon))
    print('isinstance(obj, Polygon) %s' % isinstance(obj, Polygon))


