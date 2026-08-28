from abc import ABC, abstractmethod

class _MutationalModel(ABC):
    name = ''

    @abstractmethod
    def mutate(sequence):
        return sequence

    @classmethod
    def get_subclasses(cls):
        """
        Get all subclasses of the _MutationalModel class.
        :return: A generator of all subclasses of this class.
        """
        for klass in cls.__subclasses__():
            yield klass
            for subklass in klass.get_subclasses():
                yield subklass
    @classmethod
    def find_by_name(cls, name):
        for klass in cls.get_subclasses():
            if klass.name == name:
                return klass
        return None
    
    def __repr__(self):
        return f"<{self.__class__.__name__}()>"

