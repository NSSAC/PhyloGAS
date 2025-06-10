from collections import UserDict

from .base import _MutationalModel
from .poor import *
from .simple import *
from .rate_limited import *

class _Registry(UserDict):

    def __getitem__(self, key):
        """Gets mutational model based on name"""

        # First use manually recorded models
        try:
            return super().__getitem__(key)
        except KeyError:
            klass = _MutationalModel.find_by_name(key)
            if klass is None:
                raise KeyError
            return klass
        

registry = _Registry()








