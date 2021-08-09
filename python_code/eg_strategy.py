from abc import ABC, abstractmethod

# https://refactoring.guru/design-patterns/strategy


# The strategy interface declares operations common to all
# supported versions of some algorithm. The context uses this
# interface to call the algorithm defined by the concrete
# strategies.
class Strategy(ABC):
    @abstractmethod
    def execute(self, rotor):
        pass


# Concrete strategies implement the algorithm while following
# the base strategy interface. The interface makes them
# interchangeable in the context.
class ConcreteStrategyAdd(Strategy):
    def execute(self, rotor):
        print('[ConcreteStrategyAdd.execute]')
        return rotor + rotor


class ConcreteStrategySubtract(Strategy):
    def execute(self, rotor):
        print('[ConcreteStrategySubtract.execute]')
        return None


class ConcreteStrategyMultiply(Strategy):
    def execute(self, rotor):
        print('[ConcreteStrategyMultiply.execute]')
        return rotor + '\n' + rotor


# The context defines the interface of interest to clients.
class Context:
    def __init__(self, strategy):
        # # The context maintains a reference to one of the strategy
        # # objects. The context doesn't know the concrete class of a
        # # strategy. It should work with all strategies via the
        # # strategy interface.
        # private strategy: Strategy
        #
        # # Usually the context accepts a strategy through the
        # # constructor, and also provides a setter so that the
        # # strategy can be switched at runtime.
        # method setStrategy(Strategy strategy) is
        #     this.strategy = strategy
        assert (isinstance(strategy, Strategy))
        self.strategy = strategy

    # The context delegates some work to the strategy object
    # instead of implementing multiple versions of the
    # algorithm on its own.
    def executeStrategy(self, rotor):
        return self.strategy.execute(rotor)


# The client code picks a concrete strategy and passes it to
# the context. The client should be aware of the differences
# between strategies in order to make the right choice.
class ExampleApplication:
    def __init__(self, action):
        # Create context object.
        # Read first number.
        # Read last number.
        # Read the desired action from user input.

        self.context = None
        if action == 'addition':
            self.context = Context(ConcreteStrategyAdd())
        if action == "subtraction":
            self.context = Context(ConcreteStrategySubtract())
        if action == "multiplication":
            self.context = Context(ConcreteStrategyMultiply())

        rotor = """
        \           /
         \====x====/
          |       |
        """
        result = self.context.executeStrategy(rotor)

        # Print result.
        print('[ExampleApplication] result %s' % result)


if __name__ == "__main__":
    app1 = ExampleApplication('addition')
    app2 = ExampleApplication('subtraction')
    app3 = ExampleApplication('multiplication')
